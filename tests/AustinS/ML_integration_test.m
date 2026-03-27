%% 5G NR DLSCH HARQ + Reinforcement Learning (Single Consolidated Script)
% Adaptive HARQ Redundancy Version (RV) selection using MATLAB RL Toolbox
%
% This script can:
%   1) Train a DQN agent to choose HARQ RV adaptively
%   2) Evaluate the trained agent in the 5G NR DLSCH/PDSCH link simulation
%   3) Compare against a fixed RV baseline
%
% Requirements:
%   - 5G Toolbox
%   - Reinforcement Learning Toolbox
%   - Existing helper functions/classes already used in your original script:
%       HARQEntity
%       getInitialChannelEstimate
%       getPrecodingMatrix
%       generateAWGN
%       hSkipWeakTimingOffset
%       precodeChannelEstimate
%       getConstellationRefPoints
%       evmMerMetricsDD
%       logPerfMetrics
%
% Author note:
%   This is a practical first-pass integration where RL selects the RV per slot.
%   The PHY chain remains intact. The agent action is one of RV = {0,1,2,3}.

clear; clc;

%% ========================= USER CONFIGURATION =========================
% --- Top-level mode controls ---
TRAIN_RL          = true;    % true = train DQN agent first
SAVE_AGENT        = true;    % save trained agent
LOAD_AGENT_IF_AVL = true;    % load saved agent if present and TRAIN_RL=false
USE_RL_FOR_EVAL   = true;    % true = use trained RL agent for final evaluation
RUN_FIXED_BASELINE = true;   % true = also run baseline fixed RV simulation for comparison

agentFile = "trainedHarqDqnAgent.mat";

% --- Simulation parameters ---
cfg = struct();
cfg.SNRdB = 8;
cfg.Modulation = "QPSK";           % 'QPSK','16QAM','64QAM','256QAM','1024QAM'
cfg.NHARQProcesses = 16;           % [1 2 4 8 16]
cfg.rvSeqBaseline = [0 2 3 1];     % fixed baseline sequence
cfg.nTxAnts = 8;
cfg.nRxAnts = 8;
cfg.NumLayers = 2;
cfg.DelayProfile = "TDL-C";
cfg.totalNoSlots = 20;
cfg.perfectEstimation = false;
cfg.dispPlots = 0;
cfg.logging = 1;

% --- RL training parameters ---
rlcfg = struct();
rlcfg.MaxEpisodes = 300;
rlcfg.MaxStepsPerEpisode = cfg.totalNoSlots;
rlcfg.ScoreWindow = 20;
rlcfg.StopTrainingValue = 60;   % tune as needed
rlcfg.UseDoubleDQN = true;
rlcfg.MiniBatchSize = 128;
rlcfg.ExperienceBufferLength = 1e5;
rlcfg.DiscountFactor = 0.99;
rlcfg.TargetSmoothFactor = 1e-3;
rlcfg.InitialEpsilon = 1.0;
rlcfg.EpsilonDecay = 1e-4;
rlcfg.EpsilonMin = 0.05;

rng("default");

%% ========================= OPTIONAL LOGGING SETUP =========================
if cfg.logging == 1
    logDir = fullfile(pwd, 'logs');
    if ~exist(logDir,'dir'); mkdir(logDir); end

    figDir = fullfile(pwd,'figures');
    if ~exist(figDir,'dir'); mkdir(figDir); end

    rvStr = sprintf('%d', cfg.rvSeqBaseline);
    antStr = sprintf('%dx%d', cfg.nTxAnts, cfg.nRxAnts);

    masterLogFile = fullfile(logDir, ...
        sprintf("HARQ_RL_%s_SNR%.1f_%s_NHARQ%d_rvSeq[%s]_%s_%s.txt", ...
        cfg.Modulation, cfg.SNRdB, cfg.DelayProfile, cfg.NHARQProcesses, ...
        rvStr, antStr, datestr(now,'yyyymmdd_HHMMSS')));

    diary(masterLogFile);
end

fprintf('\n=========================================================\n');
fprintf('     5G NR DLSCH HARQ + Reinforcement Learning Script\n');
fprintf('=========================================================\n');
disp(cfg);

%% ========================= TRAIN OR LOAD RL AGENT =========================
agent = [];

if TRAIN_RL
    fprintf('\n=== TRAINING RL AGENT ===\n');

    % Observation: [meanSINR, avgGain, log10(condNum), isNewTB, harqRound, prevAck, prevRV, configuredSNR]
    obsInfo = rlNumericSpec([8 1], ...
        LowerLimit=-inf*ones(8,1), ...
        UpperLimit= inf*ones(8,1));
    obsInfo.Name = "HARQState";

    % Action: choose RV in {0,1,2,3}
    actInfo = rlFiniteSetSpec({0,1,2,3});
    actInfo.Name = "RVAction";

    env = rlFunctionEnv(obsInfo, actInfo, ...
        @(action,logged) localStepFunction(action, logged), ...
        @() localResetFunction(cfg));

    % DQN network
    statePath = [
        featureInputLayer(obsInfo.Dimension(1), Name="state")
        fullyConnectedLayer(64, Name="fc1")
        reluLayer(Name="relu1")
        fullyConnectedLayer(64, Name="fc2")
        reluLayer(Name="relu2")
        fullyConnectedLayer(numel(actInfo.Elements), Name="qOut")
    ];

    criticNet = dlnetwork(layerGraph(statePath));
    critic = rlVectorQValueFunction(criticNet, obsInfo, actInfo);

    agentOpts = rlDQNAgentOptions( ...
        UseDoubleDQN=rlcfg.UseDoubleDQN, ...
        TargetSmoothFactor=rlcfg.TargetSmoothFactor, ...
        ExperienceBufferLength=rlcfg.ExperienceBufferLength, ...
        MiniBatchSize=rlcfg.MiniBatchSize, ...
        DiscountFactor=rlcfg.DiscountFactor);

    agentOpts.EpsilonGreedyExploration.Epsilon    = rlcfg.InitialEpsilon;
    agentOpts.EpsilonGreedyExploration.EpsilonDecay = rlcfg.EpsilonDecay;
    agentOpts.EpsilonGreedyExploration.EpsilonMin   = rlcfg.EpsilonMin;

    agent = rlDQNAgent(critic, agentOpts);

    trainOpts = rlTrainingOptions( ...
        MaxEpisodes=rlcfg.MaxEpisodes, ...
        MaxStepsPerEpisode=rlcfg.MaxStepsPerEpisode, ...
        ScoreAveragingWindowLength=rlcfg.ScoreWindow, ...
        Verbose=true, ...
        Plots="training-progress", ...
        StopTrainingCriteria="AverageReward", ...
        StopTrainingValue=rlcfg.StopTrainingValue);

    trainingStats = train(agent, env, trainOpts); %#ok<NASGU>

    if SAVE_AGENT
        save(agentFile, "agent", "trainingStats", "cfg", "rlcfg");
        fprintf('Saved trained agent to: %s\n', agentFile);
    end

elseif LOAD_AGENT_IF_AVL && isfile(agentFile)
    fprintf('\n=== LOADING EXISTING RL AGENT ===\n');
    S = load(agentFile, "agent");
    agent = S.agent;
    fprintf('Loaded agent from: %s\n', agentFile);
end

%% ========================= EVALUATION RUNS =========================
baselineResults = [];
rlResults = [];

if RUN_FIXED_BASELINE
    fprintf('\n=== RUNNING FIXED RV BASELINE ===\n');
    baselineResults = run5GHarqSimulation(cfg, [], false);
end

if USE_RL_FOR_EVAL
    if isempty(agent)
        warning('USE_RL_FOR_EVAL is true, but no trained/loaded agent is available. Skipping RL evaluation.');
    else
        fprintf('\n=== RUNNING RL-ADAPTIVE HARQ EVALUATION ===\n');
        rlResults = run5GHarqSimulation(cfg, agent, true);
    end
end

%% ========================= SUMMARY COMPARISON =========================
fprintf('\n=========================================================\n');
fprintf('                    FINAL COMPARISON\n');
fprintf('=========================================================\n');

if ~isempty(baselineResults)
    fprintf('\n--- Fixed RV Baseline ---\n');
    printResultsSummary(baselineResults);
end

if ~isempty(rlResults)
    fprintf('\n--- RL Adaptive RV ---\n');
    printResultsSummary(rlResults);
end

if ~isempty(baselineResults) && ~isempty(rlResults)
    fprintf('\n--- Delta (RL - Baseline) ---\n');
    fprintf('Attempt BLER delta:       %+0.4f\n', rlResults.attemptBLER - baselineResults.attemptBLER);
    fprintf('Final BLER delta:         %+0.4f\n', rlResults.finalBLER - baselineResults.finalBLER);
    fprintf('Retransmission rate delta:%+0.2f %%\n', rlResults.retxRatePct - baselineResults.retxRatePct);
    fprintf('Throughput eff delta:     %+0.2f %%\n', rlResults.throughputEfficiencyPct - baselineResults.throughputEfficiencyPct);
    fprintf('Avg throughput delta:     %+0.2f Mbps\n', rlResults.avgThroughputMbps - baselineResults.avgThroughputMbps);
    fprintf('Spectral efficiency delta:%+0.4f bits/s/Hz\n', rlResults.spectralEfficiency - baselineResults.spectralEfficiency);
end

fprintf('\nDone.\n');

if cfg.logging == 1
    diary off;
end

%% ========================================================================
%% ============================ MAIN SIM FUNCTION ==========================
%% ========================================================================
function results = run5GHarqSimulation(cfg, agent, useRL)

    % -----------------------------
    % Local run configuration
    % -----------------------------
    SNRdB = cfg.SNRdB;
    Modulation = cfg.Modulation;
    NHARQProcesses = cfg.NHARQProcesses;
    rvSeq = cfg.rvSeqBaseline;
    nTxAnts = cfg.nTxAnts;
    nRxAnts = cfg.nRxAnts;
    NumLayers = cfg.NumLayers;
    DelayProfile = cfg.DelayProfile;
    totalNoSlots = cfg.totalNoSlots;
    perfectEstimation = cfg.perfectEstimation;
    dispPlots = cfg.dispPlots;

    % -----------------------------
    % Carrier configuration
    % -----------------------------
    carrier = nrCarrierConfig;

    % -----------------------------
    % PDSCH configuration
    % -----------------------------
    pdsch = nrPDSCHConfig;
    pdsch.Modulation = Modulation;
    pdsch.NumLayers = NumLayers;
    pdsch.PRBSet = 0:carrier.NSizeGrid-1;

    pdsch.DMRS.DMRSAdditionalPosition = 1;
    pdsch.DMRS.DMRSConfigurationType = 1;
    pdsch.DMRS.DMRSLength = 2;

    if pdsch.NumCodewords == 1
        codeRate = 490/1024;
    else
        codeRate = [490 490]./1024;
    end

    % -----------------------------
    % DL-SCH encoder/decoder
    % -----------------------------
    encodeDLSCH = nrDLSCH;
    encodeDLSCH.MultipleHARQProcesses = true;
    encodeDLSCH.TargetCodeRate = codeRate;

    decodeDLSCH = nrDLSCHDecoder;
    decodeDLSCH.MultipleHARQProcesses = true;
    decodeDLSCH.TargetCodeRate = codeRate;
    decodeDLSCH.LDPCDecodingAlgorithm = "Normalized min-sum";
    decodeDLSCH.MaximumLDPCIterationCount = 6;

    % -----------------------------
    % HARQ entity
    % -----------------------------
    harqEntity = HARQEntity(0:NHARQProcesses-1, rvSeq, pdsch.NumCodewords);

    % -----------------------------
    % Channel configuration
    % -----------------------------
    if pdsch.NumLayers > min(nTxAnts,nRxAnts)
        error("The number of layers ("+string(pdsch.NumLayers)+") must be smaller than min(nTxAnts,nRxAnts) ("+string(min(nTxAnts,nRxAnts))+")");
    end

    channel = nrTDLChannel;
    channel.DelayProfile = DelayProfile;
    channel.NumTransmitAntennas = nTxAnts;
    channel.NumReceiveAntennas = nRxAnts;

    ofdmInfo = nrOFDMInfo(carrier);
    channel.SampleRate = ofdmInfo.SampleRate;
    channel.ChannelResponseOutput = 'ofdm-response';

    chInfo = info(channel);

    fprintf('\n===== RUN START =====\n');
    fprintf('Timestamp: %s\n', datestr(now));
    fprintf('Mode: %s\n', ternary(useRL, 'RL Adaptive RV', 'Fixed Baseline RV'));
    fprintf('SNR: %.2f dB | Modulation: %s | NHARQ: %d | DelayProfile: %s\n', ...
        SNRdB, Modulation, NHARQProcesses, DelayProfile);
    fprintf('Tx/Rx Antennas: %d/%d | Layers: %d | Slots: %d\n', ...
        nTxAnts, nRxAnts, NumLayers, totalNoSlots);

    fprintf('\n========== Channel Model Configuration ==========\n');
    fprintf('Channel Type: %s\n', channel.DelayProfile);
    fprintf('Sample Rate: %.2f MHz\n', channel.SampleRate / 1e6);
    fprintf('Maximum Channel Delay: %d samples\n', chInfo.MaximumChannelDelay);
    fprintf('Transmit Antennas: %d\n', channel.NumTransmitAntennas);
    fprintf('Receive Antennas: %d\n', channel.NumReceiveAntennas);

    fprintf('\n--- Path Characteristics ---\n');
    fprintf('Path Delays (ns): [%s]\n', num2str(chInfo.PathDelays * 1e9, '%.2f '));
    fprintf('Average Path Gains (dB): [%s]\n', num2str(chInfo.AveragePathGains, '%.2f '));

    pathPowers = 10.^(chInfo.AveragePathGains/10);
    meanDelay = sum(chInfo.PathDelays .* pathPowers) / sum(pathPowers);
    rmsDelaySpread = sqrt(sum(((chInfo.PathDelays - meanDelay).^2) .* pathPowers) / sum(pathPowers));

    fprintf('RMS Delay Spread: %.2f ns\n', rmsDelaySpread * 1e9);
    fprintf('Mean Excess Delay: %.2f ns\n', meanDelay * 1e9);

    if isprop(channel, 'MaximumDopplerShift')
        fprintf('Maximum Doppler Shift: %.2f Hz\n', channel.MaximumDopplerShift);
    end
    fprintf('================================================\n\n');

    % -----------------------------
    % Optional constellation diagram object
    % -----------------------------
    constPlot = comm.ConstellationDiagram;
    constPlot.ReferenceConstellation = getConstellationRefPoints(pdsch.Modulation);
    constPlot.EnableMeasurements = 1;

    % -----------------------------
    % Tracking variables
    % -----------------------------
    offset = 0;

    totalTxBits = 0;
    totalRxBits = 0;
    totalTransmissions = 0;
    totalInitialTransmissions = 0;
    totalRetransmissions = 0;
    blockErrors = 0;
    successfulBlocks = 0;

    perSlotSuccess = zeros(totalNoSlots, 1);
    perSlotBER = zeros(totalNoSlots, 1);

    estChannelGrid = getInitialChannelEstimate(channel,carrier);
    newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChannelGrid);

    instantaneousSINR = zeros(totalNoSlots, pdsch.NumLayers);
    avgChannelGain = zeros(totalNoSlots, 1);
    conditionNumber = zeros(totalNoSlots, 1);

    EVMrms_pct = zeros(totalNoSlots, pdsch.NumLayers);
    EVMpk_pct  = zeros(totalNoSlots, pdsch.NumLayers);
    EVMrms_dB  = zeros(totalNoSlots, pdsch.NumLayers);
    EVMpk_dB   = zeros(totalNoSlots, pdsch.NumLayers);
    MERavg_dB  = zeros(totalNoSlots, pdsch.NumLayers);

    slotDuration_s = 1e-3 / carrier.SlotsPerSubframe;

    perfState = struct( ...
        'totalTxBits', 0, ...
        'totalRxBits', 0, ...
        'totalTransmissions', 0, ...
        'totalInitialTransmissions', 0, ...
        'totalRetransmissions', 0, ...
        'blockErrors', 0, ...
        'successfulBlocks', 0, ...
        'attemptsTotal', 0, ...
        'attemptsFailed', 0);

    perfLog = repmat(struct( ...
        'slot', 0, ...
        'isNewTB', false, ...
        'txBitsThisSlot', 0, ...
        'rxBitsThisSlot', 0, ...
        'blkErrCountThisSlot', 0, ...
        'attemptBLER_thisSlot', 0, ...
        'attemptBLER_cum', 0, ...
        'finalBLER', 0, ...
        'totalTransmissions', 0, ...
        'totalInitialTransmissions', 0, ...
        'totalRetransmissions', 0, ...
        'successfulBlocks', 0, ...
        'blockErrors', 0, ...
        'throughputEfficiencyPct', 0, ...
        'instThroughputMbps', 0, ...
        'avgThroughputMbps', 0), totalNoSlots, 1);

    % RL feedback state
    prevAck = 1;
    prevRV = 0;
    prevMeanSINR = SNRdB;
    prevAvgChannelGain = 0;
    prevConditionNumber = 1;
    currentHarqRound = 0;

    rlSelectedRVs = zeros(totalNoSlots,1);
    rlRewards = zeros(totalNoSlots,1);

    % -----------------------------
    % MAIN LOOP
    % -----------------------------
    for nSlot = 0:totalNoSlots-1

        carrier.NSlot = nSlot;

        [pdschIndices,pdschInfo] = nrPDSCHIndices(carrier,pdsch);

        Xoh_PDSCH = 0;
        trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet), ...
            pdschInfo.NREPerPRB,codeRate,Xoh_PDSCH);

        % --- New TB handling ---
        isNewTB_any = false;
        for cwIdx = 1:pdsch.NumCodewords
            if harqEntity.NewData(cwIdx)
                isNewTB_any = true;
                trBlk = randi([0 1],trBlkSizes(cwIdx),1);
                setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);

                if harqEntity.SequenceTimeout(cwIdx)
                    resetSoftBuffer(decodeDLSCH,cwIdx-1,harqEntity.HARQProcessID);
                end
            end
        end

        if isNewTB_any
            currentHarqRound = 0;
        end

        % --- RL action or fixed RV selection ---
        if useRL && ~isempty(agent)
            obs = double([ ...
                prevMeanSINR; ...
                prevAvgChannelGain; ...
                log10(max(prevConditionNumber,1)); ...
                double(currentHarqRound == 0); ...
                currentHarqRound; ...
                prevAck; ...
                prevRV; ...
                SNRdB ]);

            rvAction = getAction(agent, {obs});
            selectedRV = rvAction{1};
        else
            selectedRV = rvSeq(mod(currentHarqRound, numel(rvSeq)) + 1);
        end

        rlSelectedRVs(nSlot+1) = selectedRV;

        rvVec = repmat(selectedRV, 1, pdsch.NumCodewords);

        codedTrBlock = encodeDLSCH(pdsch.Modulation,pdsch.NumLayers,pdschInfo.G, ...
            rvVec,harqEntity.HARQProcessID);

        pdschSymbols = nrPDSCH(carrier,pdsch,codedTrBlock);

        precodingWeights = newPrecodingWeight;
        pdschSymbolsPrecoded = pdschSymbols*precodingWeights;

        dmrsSymbols = nrPDSCHDMRS(carrier,pdsch);
        dmrsIndices = nrPDSCHDMRSIndices(carrier,pdsch);

        pdschGrid = nrResourceGrid(carrier,nTxAnts);

        [~,pdschAntIndices] = nrExtractResources(pdschIndices,pdschGrid);
        pdschGrid(pdschAntIndices) = pdschSymbolsPrecoded;

        for p = 1:size(dmrsSymbols,2)
            [~,dmrsAntIndices] = nrExtractResources(dmrsIndices(:,p),pdschGrid);
            pdschGrid(dmrsAntIndices) = pdschGrid(dmrsAntIndices) + ...
                dmrsSymbols(:,p)*precodingWeights(p,:);
        end

        [txWaveform,waveformInfo] = nrOFDMModulate(carrier,pdschGrid);

        chInfo = info(channel);
        maxChDelay = chInfo.MaximumChannelDelay;
        txWaveform = [txWaveform; zeros(maxChDelay,size(txWaveform,2))];

        [rxWaveform,ofdmChannelResponse,timingOffset] = channel(txWaveform,carrier);
        [noise,nVar] = generateAWGN(SNRdB,nRxAnts,waveformInfo.Nfft,size(rxWaveform));
        rxWaveform = rxWaveform + noise;

        if perfectEstimation
            offset = timingOffset;
        else
            [t,mag] = nrTimingEstimate(carrier,rxWaveform,dmrsIndices,dmrsSymbols);
            offset = hSkipWeakTimingOffset(offset,t,mag);
        end

        rxWaveform = rxWaveform(1+offset:end,:);
        rxGrid = nrOFDMDemodulate(carrier,rxWaveform);

        if perfectEstimation
            estChGridAnts = ofdmChannelResponse;
            noiseEst = nVar;
            newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChGridAnts);
            estChGridLayers = precodeChannelEstimate(estChGridAnts,precodingWeights.');
        else
            [estChGridLayers,noiseEst] = nrChannelEstimate(carrier,rxGrid,dmrsIndices, ...
                dmrsSymbols,'CDMLengths',pdsch.DMRS.CDMLengths);

            estChGridAnts = precodeChannelEstimate(estChGridLayers,conj(precodingWeights));
            newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChGridAnts);
        end

        % --- Channel stats ---
        channelMagnitude = abs(estChGridLayers);
        avgChannelGain(nSlot+1) = 10*log10(mean(channelMagnitude(:).^2));

        try
            H_layer = squeeze(estChGridLayers(:,:,1,:));
            H_matrix = squeeze(H_layer(round(size(H_layer,1)/2), :, :));
            if ~isempty(H_matrix)
                conditionNumber(nSlot+1) = cond(H_matrix);
            else
                conditionNumber(nSlot+1) = 1;
            end
        catch
            conditionNumber(nSlot+1) = 1;
        end

        [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxGrid,estChGridLayers);
        [pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);

        % --- EVM/MER stats ---
        refConst_col = getConstellationRefPoints(pdsch.Modulation);
        refConst_col = refConst_col(:);

        for layerIdx = 1:pdsch.NumLayers
            layerSyms = pdschEq(:, layerIdx);

            M = evmMerMetricsDD(layerSyms, refConst_col);

            EVMrms_pct(nSlot+1, layerIdx) = M.RmsEVM_pct;
            EVMpk_pct (nSlot+1, layerIdx) = M.PeakEVM_pct;
            EVMrms_dB (nSlot+1, layerIdx) = M.AvgEVM_dB;
            EVMpk_dB  (nSlot+1, layerIdx) = M.PeakEVM_dB;
            MERavg_dB (nSlot+1, layerIdx) = M.AvgMER_dB;
        end

        % --- SINR estimate per layer ---
        for layerIdx = 1:pdsch.NumLayers
            layerSymbols = pdschEq(:, layerIdx);
            [~, symbolDecisions] = min(abs(layerSymbols - refConst_col.'), [], 2);
            decidedSymbols = refConst_col(symbolDecisions);

            signalPower = mean(abs(decidedSymbols).^2);
            noisePower  = mean(abs(layerSymbols - decidedSymbols).^2);

            instantaneousSINR(nSlot+1, layerIdx) = 10*log10(max(signalPower / max(noisePower,1e-12), 1e-12));
        end

        if dispPlots
            mesh(abs(estChGridLayers(:,:,1,1)));
            title('Channel Estimate');
            xlabel('OFDM Symbol');
            ylabel("Subcarrier");
            zlabel("Magnitude");

            constPlot.ChannelNames = "Layer "+(pdsch.NumLayers:-1:1);
            constPlot.ShowLegend = true;
            constPlot(fliplr(pdschEq));
        end

        [dlschLLRs,rxSymbols] = nrPDSCHDecode(carrier,pdsch,pdschEq,noiseEst);

        csi = nrLayerDemap(csi);
        for cwIdx = 1:pdsch.NumCodewords
            Qm = length(dlschLLRs{cwIdx})/length(rxSymbols{cwIdx});
            csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);
            dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:);
        end

        decodeDLSCH.TransportBlockLength = trBlkSizes;
        [~,blkerr] = decodeDLSCH(dlschLLRs,pdsch.Modulation,pdsch.NumLayers, ...
            rvVec,harqEntity.HARQProcessID);

        statusReport = updateAndAdvance(harqEntity,blkerr,trBlkSizes,pdschInfo.G); %#ok<NASGU>

        % --- Per-slot logging ---
        [perfState, perfLog(nSlot+1)] = logPerfMetrics(perfState, nSlot, trBlkSizes, blkerr, harqEntity, slotDuration_s);

        % --- Global metrics ---
        totalTransmissions = totalTransmissions + 1;

        for cwIdx = 1:pdsch.NumCodewords
            if harqEntity.NewData(cwIdx)
                totalInitialTransmissions = totalInitialTransmissions + 1;
                totalTxBits = totalTxBits + trBlkSizes(cwIdx);
            else
                totalRetransmissions = totalRetransmissions + 1;
            end

            if ~blkerr(cwIdx)
                successfulBlocks = successfulBlocks + 1;
                totalRxBits = totalRxBits + trBlkSizes(cwIdx);
                perSlotSuccess(nSlot+1) = 1;
            else
                if harqEntity.SequenceTimeout(cwIdx)
                    blockErrors = blockErrors + 1;
                end
            end
        end

        perSlotBER(nSlot+1) = sum(blkerr) / pdsch.NumCodewords;

        % --- RL-style reward tracking even during eval ---
        ackSuccess = all(~blkerr);
        finalFailure = any(blkerr) && any(harqEntity.SequenceTimeout);
        txBitsThisSlot = sum(trBlkSizes);
        rxBitsThisSlot = double(ackSuccess) * txBitsThisSlot;
        isRetransmission = ~(currentHarqRound == 0);

        reward = ...
            5.0 * double(ackSuccess) ...
            - 3.0 * double(finalFailure) ...
            - 0.25 * double(isRetransmission) ...
            + 1e-4 * double(rxBitsThisSlot) ...
            - 5e-5 * double(txBitsThisSlot);

        rlRewards(nSlot+1) = reward;

        fprintf(['Slot %3d | Mode=%s | RV=%d | HarqRound=%d | NewTB=%d | ' ...
                 'AttemptBLER=%.3f | CumAttemptBLER=%.3f | FinalBLER=%.3f | ' ...
                 'Eff=%.1f%% | InstThr=%.2f | AvgThr=%.2f | Reward=%.3f\n'], ...
            nSlot, ternary(useRL,'RL','FIXED'), selectedRV, currentHarqRound, ...
            perfLog(nSlot+1).isNewTB, ...
            perfLog(nSlot+1).attemptBLER_thisSlot, perfLog(nSlot+1).attemptBLER_cum, ...
            perfLog(nSlot+1).finalBLER, perfLog(nSlot+1).throughputEfficiencyPct, ...
            perfLog(nSlot+1).instThroughputMbps, perfLog(nSlot+1).avgThroughputMbps, reward);

        % --- Feedback for next slot ---
        prevAck = double(ackSuccess);
        prevRV = selectedRV;
        prevMeanSINR = mean(instantaneousSINR(nSlot+1,:));
        prevAvgChannelGain = avgChannelGain(nSlot+1);
        prevConditionNumber = max(conditionNumber(nSlot+1),1);

        if ackSuccess
            currentHarqRound = 0;
        else
            currentHarqRound = currentHarqRound + 1;
        end

        lastEstChGridLayers = estChGridLayers; %#ok<NASGU>
        lastPdschEq = pdschEq; %#ok<NASGU>
    end

    % --- Final metrics ---
    finalBLER = blockErrors / max(1,totalInitialTransmissions);
    attemptBLER = perfState.attemptsFailed / max(1, perfState.attemptsTotal);

    avgBitsPerSlot = totalRxBits / totalNoSlots;

    slotDuration = carrier.SlotsPerSubframe / (carrier.SubcarrierSpacing / 15e3) * 1e-3;
    throughput_Mbps = (avgBitsPerSlot / slotDuration) / 1e6;

    numRBs = numel(pdsch.PRBSet);
    bandwidth_Hz = numRBs * 12 * carrier.SubcarrierSpacing * 1e3;
    spectralEfficiency = (avgBitsPerSlot / slotDuration) / bandwidth_Hz;

    retxRatePct = (totalRetransmissions / max(1,totalTransmissions)) * 100;
    throughputEfficiencyPct = (totalRxBits / max(1,totalTxBits)) * 100;

    fprintf('\n========================================\n');
    fprintf('     5G NR HARQ Simulation Results      \n');
    fprintf('========================================\n');

    fprintf('\n--- Configuration ---\n');
    fprintf('HARQ Type: Chase Combining\n');
    fprintf('Mode: %s\n', ternary(useRL, 'RL Adaptive RV', 'Fixed RV Sequence'));
    fprintf('Baseline RV Sequence: [%s]\n', num2str(rvSeq));
    fprintf('Number of HARQ Processes: %d\n', NHARQProcesses);
    fprintf('Modulation: %s\n', pdsch.Modulation);
    fprintf('Number of Layers: %d\n', pdsch.NumLayers);
    fprintf('Code Rate: %.4f\n', codeRate);
    fprintf('SNR: %.2f dB\n', SNRdB);
    fprintf('Number of Slots: %d\n', totalNoSlots);

    fprintf('\n--- Transmission Statistics ---\n');
    fprintf('Total Transmissions: %d\n', totalTransmissions);
    fprintf('Initial Transmissions: %d\n', totalInitialTransmissions);
    fprintf('Retransmissions: %d\n', totalRetransmissions);
    fprintf('Average Transmissions per TB: %.2f\n', totalTransmissions / max(1,totalInitialTransmissions));
    fprintf('Retransmission Rate: %.2f%%\n', retxRatePct);

    fprintf('\n--- Performance Metrics ---\n');
    fprintf('Successful Blocks: %d / %d\n', successfulBlocks, totalInitialTransmissions);
    fprintf('Failed Blocks (after max retx): %d\n', blockErrors);
    fprintf('Attempt BLER (pre-HARQ): %.4f (%.2f%%)\n', attemptBLER, attemptBLER*100);
    fprintf('Final BLER (post-HARQ):  %.4f (%.2f%%)\n', finalBLER, finalBLER*100);

    fprintf('\n--- Throughput Analysis ---\n');
    fprintf('Total Bits Transmitted: %d\n', totalTxBits);
    fprintf('Total Bits Received (success): %d\n', totalRxBits);
    fprintf('Throughput Efficiency: %.2f%%\n', throughputEfficiencyPct);
    fprintf('Effective Code Rate: %.4f\n', totalRxBits / max(1,totalTransmissions * pdschInfo.G));
    fprintf('Average Throughput: %.2f Mbps\n', throughput_Mbps);
    fprintf('Average Bits per Slot: %.0f bits\n', avgBitsPerSlot);

    fprintf('\n--- Spectral Efficiency ---\n');
    fprintf('Spectral Efficiency: %.4f bits/s/Hz\n', spectralEfficiency);

    fprintf('\n--- Channel Model Details ---\n');
    fprintf('Delay Profile: %s\n', channel.DelayProfile);
    fprintf('RMS Delay Spread: %.2f ns\n', rmsDelaySpread * 1e9);
    fprintf('Mean Excess Delay: %.2f ns\n', meanDelay * 1e9);
    fprintf('Maximum Channel Delay: %d samples (%.2f us)\n', ...
        chInfo.MaximumChannelDelay, chInfo.MaximumChannelDelay / channel.SampleRate * 1e6);

    fprintf('\n--- Channel Quality Statistics ---\n');
    fprintf('Average Channel Gain: %.2f dB\n', mean(avgChannelGain));
    fprintf('Channel Gain Std Dev: %.2f dB\n', std(avgChannelGain));
    fprintf('Average Condition Number: %.2f\n', mean(conditionNumber(conditionNumber > 0)));

    fprintf('\n--- SINR Statistics (per layer) ---\n');
    for layerIdx = 1:pdsch.NumLayers
        fprintf('Layer %d - Mean SINR: %.2f dB, Std: %.2f dB\n', ...
            layerIdx, mean(instantaneousSINR(:, layerIdx)), std(instantaneousSINR(:, layerIdx)));
    end

    coherenceBW_Hz = 1 / (5 * rmsDelaySpread);
    fprintf('\nCoherence Bandwidth (50%% corr): %.2f kHz\n', coherenceBW_Hz / 1e3);

    if isprop(channel, 'MaximumDopplerShift') && channel.MaximumDopplerShift > 0
        coherenceTime_s = 9 / (16 * pi * channel.MaximumDopplerShift);
        fprintf('Coherence Time (50%% corr): %.2f ms\n', coherenceTime_s * 1e3);
    end

    if useRL
        fprintf('\n--- RL Action Statistics ---\n');
        fprintf('Mean selected RV: %.2f\n', mean(rlSelectedRVs));
        fprintf('RV0 count: %d\n', sum(rlSelectedRVs==0));
        fprintf('RV1 count: %d\n', sum(rlSelectedRVs==1));
        fprintf('RV2 count: %d\n', sum(rlSelectedRVs==2));
        fprintf('RV3 count: %d\n', sum(rlSelectedRVs==3));
        fprintf('Mean reward per slot: %.3f\n', mean(rlRewards));
    end

    results = struct();
    results.mode = ternary(useRL, 'RL', 'Fixed');
    results.totalTransmissions = totalTransmissions;
    results.totalInitialTransmissions = totalInitialTransmissions;
    results.totalRetransmissions = totalRetransmissions;
    results.successfulBlocks = successfulBlocks;
    results.blockErrors = blockErrors;
    results.attemptBLER = attemptBLER;
    results.finalBLER = finalBLER;
    results.totalTxBits = totalTxBits;
    results.totalRxBits = totalRxBits;
    results.throughputEfficiencyPct = throughputEfficiencyPct;
    results.avgThroughputMbps = throughput_Mbps;
    results.spectralEfficiency = spectralEfficiency;
    results.retxRatePct = retxRatePct;
    results.avgChannelGain_dB = mean(avgChannelGain);
    results.avgConditionNumber = mean(conditionNumber(conditionNumber > 0));
    results.meanSINR_dB = mean(instantaneousSINR(:));
    results.rvSelections = rlSelectedRVs;
    results.rewards = rlRewards;
end

%% ========================================================================
%% ============================ RL RESET FUNCTION ==========================
%% ========================================================================
function [initialObs, logged] = localResetFunction(cfg)

    logged = struct();

    logged.cfg = cfg;
    logged.SNRdB = cfg.SNRdB;
    logged.Modulation = cfg.Modulation;
    logged.NHARQProcesses = cfg.NHARQProcesses;
    logged.nTxAnts = cfg.nTxAnts;
    logged.nRxAnts = cfg.nRxAnts;
    logged.NumLayers = cfg.NumLayers;
    logged.DelayProfile = cfg.DelayProfile;
    logged.totalNoSlots = cfg.totalNoSlots;
    logged.nSlot = 0;
    logged.perfectEstimation = cfg.perfectEstimation;

    carrier = nrCarrierConfig;

    pdsch = nrPDSCHConfig;
    pdsch.Modulation = logged.Modulation;
    pdsch.NumLayers = logged.NumLayers;
    pdsch.PRBSet = 0:carrier.NSizeGrid-1;
    pdsch.DMRS.DMRSAdditionalPosition = 1;
    pdsch.DMRS.DMRSConfigurationType = 1;
    pdsch.DMRS.DMRSLength = 2;

    if pdsch.NumCodewords == 1
        codeRate = 490/1024;
    else
        codeRate = [490 490]./1024;
    end

    encodeDLSCH = nrDLSCH;
    encodeDLSCH.MultipleHARQProcesses = true;
    encodeDLSCH.TargetCodeRate = codeRate;

    decodeDLSCH = nrDLSCHDecoder;
    decodeDLSCH.MultipleHARQProcesses = true;
    decodeDLSCH.TargetCodeRate = codeRate;
    decodeDLSCH.LDPCDecodingAlgorithm = "Normalized min-sum";
    decodeDLSCH.MaximumLDPCIterationCount = 6;

    baseRvSeq = [0 2 3 1];
    harqEntity = HARQEntity(0:logged.NHARQProcesses-1, baseRvSeq, pdsch.NumCodewords);

    channel = nrTDLChannel;
    channel.DelayProfile = logged.DelayProfile;
    channel.NumTransmitAntennas = logged.nTxAnts;
    channel.NumReceiveAntennas = logged.nRxAnts;

    ofdmInfo = nrOFDMInfo(carrier);
    channel.SampleRate = ofdmInfo.SampleRate;
    channel.ChannelResponseOutput = "ofdm-response";

    estChannelGrid = getInitialChannelEstimate(channel,carrier);
    newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChannelGrid);

    logged.carrier = carrier;
    logged.pdsch = pdsch;
    logged.codeRate = codeRate;
    logged.encodeDLSCH = encodeDLSCH;
    logged.decodeDLSCH = decodeDLSCH;
    logged.harqEntity = harqEntity;
    logged.channel = channel;
    logged.newPrecodingWeight = newPrecodingWeight;
    logged.offset = 0;

    logged.prevAck = 1;
    logged.prevRV = 0;
    logged.prevSINR = logged.SNRdB;
    logged.prevGain = 0;
    logged.prevCond = 1;
    logged.harqRound = 0;

    logged.totalTxBits = 0;
    logged.totalRxBits = 0;
    logged.totalRetransmissions = 0;
    logged.blockErrors = 0;
    logged.successfulBlocks = 0;

    initialObs = buildObservation(logged);
end

%% ========================================================================
%% ============================ RL STEP FUNCTION ===========================
%% ========================================================================
function [nextObs,reward,isDone,logged] = localStepFunction(action, logged)

    selectedRV = action{1};

    carrier = logged.carrier;
    pdsch = logged.pdsch;
    encodeDLSCH = logged.encodeDLSCH;
    decodeDLSCH = logged.decodeDLSCH;
    harqEntity = logged.harqEntity;
    channel = logged.channel;
    newPrecodingWeight = logged.newPrecodingWeight;
    offset = logged.offset;
    SNRdB = logged.SNRdB;

    nSlot = logged.nSlot;
    carrier.NSlot = nSlot;

    [pdschIndices,pdschInfo] = nrPDSCHIndices(carrier,pdsch);
    Xoh_PDSCH = 0;
    trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet), ...
        pdschInfo.NREPerPRB,logged.codeRate,Xoh_PDSCH);

    isNewTB_any = false;
    for cwIdx = 1:pdsch.NumCodewords
        if harqEntity.NewData(cwIdx)
            isNewTB_any = true;
            trBlk = randi([0 1],trBlkSizes(cwIdx),1);
            setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);

            if harqEntity.SequenceTimeout(cwIdx)
                resetSoftBuffer(decodeDLSCH,cwIdx-1,harqEntity.HARQProcessID);
            end
        end
    end

    if isNewTB_any
        logged.harqRound = 0;
    else
        logged.harqRound = logged.harqRound + 1;
        logged.totalRetransmissions = logged.totalRetransmissions + 1;
    end

    rvVec = repmat(selectedRV,1,pdsch.NumCodewords);

    codedTrBlock = encodeDLSCH(pdsch.Modulation,pdsch.NumLayers,pdschInfo.G, ...
        rvVec,harqEntity.HARQProcessID);

    pdschSymbols = nrPDSCH(carrier,pdsch,codedTrBlock);
    precodingWeights = newPrecodingWeight;
    pdschSymbolsPrecoded = pdschSymbols*precodingWeights;

    dmrsSymbols = nrPDSCHDMRS(carrier,pdsch);
    dmrsIndices = nrPDSCHDMRSIndices(carrier,pdsch);

    pdschGrid = nrResourceGrid(carrier,logged.nTxAnts);
    [~,pdschAntIndices] = nrExtractResources(pdschIndices,pdschGrid);
    pdschGrid(pdschAntIndices) = pdschSymbolsPrecoded;

    for p = 1:size(dmrsSymbols,2)
        [~,dmrsAntIndices] = nrExtractResources(dmrsIndices(:,p),pdschGrid);
        pdschGrid(dmrsAntIndices) = pdschGrid(dmrsAntIndices) + ...
            dmrsSymbols(:,p)*precodingWeights(p,:);
    end

    [txWaveform,waveformInfo] = nrOFDMModulate(carrier,pdschGrid);

    chInfo = info(channel);
    maxChDelay = chInfo.MaximumChannelDelay;
    txWaveform = [txWaveform; zeros(maxChDelay,size(txWaveform,2))];

    [rxWaveform,ofdmChannelResponse,timingOffset] = channel(txWaveform,carrier);
    [noise,nVar] = generateAWGN(SNRdB,logged.nRxAnts,waveformInfo.Nfft,size(rxWaveform));
    rxWaveform = rxWaveform + noise;

    if logged.perfectEstimation
        offset = timingOffset;
    else
        [t,mag] = nrTimingEstimate(carrier,rxWaveform,dmrsIndices,dmrsSymbols);
        offset = hSkipWeakTimingOffset(offset,t,mag);
    end

    rxWaveform = rxWaveform(1+offset:end,:);
    rxGrid = nrOFDMDemodulate(carrier,rxWaveform);

    if logged.perfectEstimation
        estChGridAnts = ofdmChannelResponse;
        noiseEst = nVar;
        newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChGridAnts);
        estChGridLayers = precodeChannelEstimate(estChGridAnts,precodingWeights.');
    else
        [estChGridLayers,noiseEst] = nrChannelEstimate(carrier,rxGrid,dmrsIndices, ...
            dmrsSymbols,'CDMLengths',pdsch.DMRS.CDMLengths);
        estChGridAnts = precodeChannelEstimate(estChGridLayers,conj(precodingWeights));
        newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChGridAnts);
    end

    [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxGrid,estChGridLayers);
    [pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);

    [dlschLLRs,rxSymbols] = nrPDSCHDecode(carrier,pdsch,pdschEq,noiseEst);

    csi = nrLayerDemap(csi);
    for cwIdx = 1:pdsch.NumCodewords
        Qm = length(dlschLLRs{cwIdx})/length(rxSymbols{cwIdx});
        csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);
        dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:);
    end

    decodeDLSCH.TransportBlockLength = trBlkSizes;
    [~,blkerr] = decodeDLSCH(dlschLLRs,pdsch.Modulation,pdsch.NumLayers, ...
        rvVec,harqEntity.HARQProcessID);

    updateAndAdvance(harqEntity,blkerr,trBlkSizes,pdschInfo.G);

    channelMagnitude = abs(estChGridLayers);
    avgGain = 10*log10(mean(channelMagnitude(:).^2));

    condNum = 1;
    try
        H_layer = squeeze(estChGridLayers(:,:,1,:));
        H_matrix = squeeze(H_layer(round(size(H_layer,1)/2), :, :));
        if ~isempty(H_matrix)
            condNum = cond(H_matrix);
        end
    catch
        condNum = 1;
    end

    refConst = getConstellationRefPoints(pdsch.Modulation);
    refConst = refConst(:);

    sinrPerLayer = zeros(pdsch.NumLayers,1);
    for layerIdx = 1:pdsch.NumLayers
        layerSymbols = pdschEq(:,layerIdx);
        [~,idx] = min(abs(layerSymbols - refConst.'), [], 2);
        decidedSymbols = refConst(idx);
        signalPower = mean(abs(decidedSymbols).^2);
        noisePower = mean(abs(layerSymbols - decidedSymbols).^2);
        sinrPerLayer(layerIdx) = 10*log10(max(signalPower/max(noisePower,1e-12),1e-12));
    end

    meanSINR = mean(sinrPerLayer);

    ackSuccess = all(~blkerr);
    finalFailure = any(blkerr) && any(harqEntity.SequenceTimeout);

    txBitsThisSlot = sum(trBlkSizes);
    rxBitsThisSlot = double(ackSuccess) * txBitsThisSlot;
    isRetransmission = ~(logged.harqRound == 0);

    reward = ...
        5.0 * double(ackSuccess) ...
      - 3.0 * double(finalFailure) ...
      - 0.25 * double(isRetransmission) ...
      + 1e-4 * double(rxBitsThisSlot) ...
      - 5e-5 * double(txBitsThisSlot);

    logged.prevAck = double(ackSuccess);
    logged.prevRV = selectedRV;
    logged.prevSINR = meanSINR;
    logged.prevGain = avgGain;
    logged.prevCond = condNum;
    logged.newPrecodingWeight = newPrecodingWeight;
    logged.offset = offset;
    logged.harqEntity = harqEntity;
    logged.nSlot = logged.nSlot + 1;

    logged.totalTxBits = logged.totalTxBits + txBitsThisSlot;
    logged.totalRxBits = logged.totalRxBits + rxBitsThisSlot;
    logged.successfulBlocks = logged.successfulBlocks + double(ackSuccess);
    logged.blockErrors = logged.blockErrors + double(finalFailure);

    if ackSuccess
        logged.harqRound = 0;
    end

    nextObs = buildObservation(logged);
    isDone = logged.nSlot >= logged.totalNoSlots;
end

%% ========================================================================
%% ============================= OBSERVATION ===============================
%% ========================================================================
function obs = buildObservation(logged)

    isNewTB = double(logged.harqRound == 0);

    obs = [
        logged.prevSINR
        logged.prevGain
        log10(max(logged.prevCond,1))
        isNewTB
        logged.harqRound
        logged.prevAck
        logged.prevRV
        logged.SNRdB
    ];

    obs = double(obs(:));
end

%% ========================================================================
%% ============================= SUMMARY PRINT =============================
%% ========================================================================
function printResultsSummary(R)
    fprintf('Mode: %s\n', R.mode);
    fprintf('Attempt BLER: %.4f\n', R.attemptBLER);
    fprintf('Final BLER: %.4f\n', R.finalBLER);
    fprintf('Retransmission Rate: %.2f %%\n', R.retxRatePct);
    fprintf('Throughput Efficiency: %.2f %%\n', R.throughputEfficiencyPct);
    fprintf('Average Throughput: %.2f Mbps\n', R.avgThroughputMbps);
    fprintf('Spectral Efficiency: %.4f bits/s/Hz\n', R.spectralEfficiency);
    fprintf('Mean Channel Gain: %.2f dB\n', R.avgChannelGain_dB);
    fprintf('Mean SINR: %.2f dB\n', R.meanSINR_dB);
end

%% ========================================================================
%% =============================== UTILITIES ===============================
%% ========================================================================
function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end