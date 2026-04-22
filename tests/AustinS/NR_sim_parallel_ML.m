
function result = NR_sim_parallel_ML(cfg)
% NR_sim_parallel_ML(cfg)
% ML-integrated version of NR_sim_parallel
%
% Added cfg fields:
%   controlMode      : 'fixed' or 'ml'
%   fixedModulation  : '16QAM' or '64QAM' (used when controlMode='fixed')
%   mlModel          : struct with fields beta, mu, sigma
%   mlRiskThreshold  : scalar in [0,1], e.g. 0.35
%   datasetMode      : true/false, save TB-level dataset for training
%   totalNoSlots     : optional, defaults to 2000
%
% Notes:
% - Modulation is selected only on NEW transport blocks.
% - Retransmissions keep the same modulation for the TB/HARQ process.
% - Dataset rows are collected from first attempts only:
%       features = [SNRdB, Is64QAM, RecentNACKRate]
%       label    = NeedRetx (1 if blkerr on first attempt, else 0)

%% ------------------------------------------------------------------------
% Unpack config
% -------------------------------------------------------------------------
SNRdB          = cfg.SNRdB;
NHARQProcesses = cfg.NHARQProcesses;
rvSeq          = cfg.rvSeq;
nTxAnts        = cfg.nTxAnts;
nRxAnts        = cfg.nRxAnts;
NumLayers      = cfg.NumLayers;
DelayProfile   = cfg.DelayProfile;
dispPlots      = cfg.dispPlots;
logging        = cfg.logging;

if isfield(cfg,'controlMode');      controlMode = string(cfg.controlMode); else; controlMode = "fixed"; end
if isfield(cfg,'fixedModulation');  fixedModulation = string(cfg.fixedModulation); else; fixedModulation = "16QAM"; end
if isfield(cfg,'datasetMode');      datasetMode = logical(cfg.datasetMode); else; datasetMode = false; end
if isfield(cfg,'mlRiskThreshold');  mlRiskThreshold = cfg.mlRiskThreshold; else; mlRiskThreshold = 0.35; end
if isfield(cfg,'mlModel');          mlModel = cfg.mlModel; else; mlModel = []; end
if isfield(cfg,'totalNoSlots');     totalNoSlots = cfg.totalNoSlots; else; totalNoSlots = 2000; end

rng(cfg.seed, 'twister');
runTimer = tic;

%% ------------------------------------------------------------------------
% Result struct + safe defaults
% -------------------------------------------------------------------------
result = struct();
result.runIdx = cfg.runIdx;
result.logFile = "";
result.datasetFile = "";
result.ok = false;
result.skipped = false;
result.errorMessage = "";
result.controlMode = controlMode;
result.SNRdB = SNRdB;
result.NHARQProcesses = NHARQProcesses;
result.DelayProfile = string(DelayProfile);
result.ThroughputEfficiencyPct = nan;
result.HARQEfficiencyPct = nan;
result.RetransmissionRatePct = nan;
result.AttemptBLERPct = nan;
result.FinalBLERPct = nan;
result.AverageThroughputMbps = nan;
result.EstimatedLatencyMs = nan;
result.Modulation = "";

fid = -1;
logFile = "";
datasetFile = "";
logBuffer = strings(0,1);
figFileConst = "";

%% ------------------------------------------------------------------------
% Logging / output setup
% -------------------------------------------------------------------------
rvStr  = sprintf('%d', rvSeq);
antStr = sprintf('%dx%d', nTxAnts, nRxAnts);
ts = char(datetime('now','Format','yyyyMMdd_HHmmss'));

if controlMode == "fixed"
    modTag = char(fixedModulation);
else
    modTag = 'ML_16QAM_64QAM';
end

baseName = sprintf('%s_SNR%.1f_%s_NHARQ%d_rvSeq[%s]_%s_%s_run%05d', ...
    modTag, SNRdB, DelayProfile, NHARQProcesses, rvStr, antStr, ts, cfg.runIdx);

dsDir = fullfile(pwd,'ml_dataset');
if datasetMode && ~exist(dsDir,'dir'); mkdir(dsDir); end
datasetFile = fullfile(dsDir, baseName + "_dataset.mat");
result.datasetFile = string(datasetFile);

if logging == 1
    logDir = fullfile(pwd, 'logs');
    if ~exist(logDir,'dir'); mkdir(logDir); end

    figDir = fullfile(pwd,'figures');
    if ~exist(figDir,'dir'); mkdir(figDir); end

    logFile = fullfile(logDir, baseName + ".txt");

    fid = fopen(logFile, 'w');
    if fid < 0
        error('Could not open log file: %s', logFile);
    end

    figFileConst = fullfile(figDir, baseName + "_Constellation.png");
    result.logFile = string(logFile);
end

try
    %% --------------------------------------------------------------------
    % Carrier configuration
    % ---------------------------------------------------------------------
    perfectEstimation = false;

    carrier = nrCarrierConfig;
    carrier.SubcarrierSpacing = 15;   % kHz

    %% --------------------------------------------------------------------
    % PDSCH and DM-RS config
    % ---------------------------------------------------------------------
    pdsch = nrPDSCHConfig;
    pdsch.Modulation = char(fixedModulation);    % overwritten per slot if ML
    pdsch.NumLayers = NumLayers;
    pdsch.PRBSet = 0:carrier.NSizeGrid-1;

    pdsch.DMRS.DMRSAdditionalPosition = 1;
    pdsch.DMRS.DMRSConfigurationType = 1;
    pdsch.DMRS.DMRSLength = 2;

    %% --------------------------------------------------------------------
    % DL-SCH configuration
    % ---------------------------------------------------------------------
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

    harqEntity = HARQEntity(0:NHARQProcesses-1,rvSeq,pdsch.NumCodewords);

    %% --------------------------------------------------------------------
    % Channel configuration
    % ---------------------------------------------------------------------
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

    %% --------------------------------------------------------------------
    % Log run parameters
    % ---------------------------------------------------------------------
    if logging == 1
        logBuffer(end+1) = sprintf('===== RUN START =====');
        logBuffer(end+1) = sprintf('Timestamp: %s', datestr(now));
        logBuffer(end+1) = "===== RUN PARAMETERS =====";
        logBuffer(end+1) = sprintf('%15s: %g', 'SNRdB', SNRdB);
        if controlMode == "fixed"
            logBuffer(end+1) = sprintf('%15s: "%s"', 'Modulation', fixedModulation);
        else
            logBuffer(end+1) = sprintf('%15s: "%s"', 'Modulation', 'ML_16QAM_64QAM');
        end
        logBuffer(end+1) = sprintf('%15s: %d', 'NHARQProcesses', NHARQProcesses);
        logBuffer(end+1) = sprintf('%15s: [%s]', 'rvSeq', strjoin(string(rvSeq), ' '));
        logBuffer(end+1) = sprintf('%15s: %d', 'nTxAnts', nTxAnts);
        logBuffer(end+1) = sprintf('%15s: %d', 'nRxAnts', nRxAnts);
        logBuffer(end+1) = sprintf('%15s: %d', 'NumLayers', NumLayers);
        logBuffer(end+1) = sprintf('%15s: "%s"', 'DelayProfile', DelayProfile);
        logBuffer(end+1) = "===========================";
        logBuffer(end+1) = "";

        logBuffer(end+1) = sprintf('%15s: "%s"', 'ControlMode', controlMode);
        if controlMode == "fixed"
            logBuffer(end+1) = sprintf('%15s: "%s"', 'FixedMod', fixedModulation);
        else
            logBuffer(end+1) = sprintf('%15s: %.3f', 'MLRiskThr', mlRiskThreshold);
        end
        logBuffer(end+1) = sprintf('%15s: %d', 'TotalSlots', totalNoSlots);
        logBuffer(end+1) = sprintf('%15s: %d', 'DatasetMode', datasetMode);
        logBuffer(end+1) = "";

        logBuffer(end+1) = sprintf('\n========== Channel Model Configuration ==========');
        logBuffer(end+1) = sprintf('Channel Type: %s', channel.DelayProfile);
        logBuffer(end+1) = sprintf('Sample Rate: %.2f MHz', channel.SampleRate / 1e6);
        logBuffer(end+1) = sprintf('Maximum Channel Delay: %d samples', chInfo.MaximumChannelDelay);
        logBuffer(end+1) = sprintf('Transmit Antennas: %d', channel.NumTransmitAntennas);
        logBuffer(end+1) = sprintf('Receive Antennas: %d', channel.NumReceiveAntennas);
        logBuffer(end+1) = sprintf('--- Path Characteristics ---\n');
        logBuffer(end+1) = sprintf('Path Delays (ns): [%s]', num2str(chInfo.PathDelays * 1e9, '%.2f '));
        logBuffer(end+1) = sprintf('Average Path Gains (dB): [%s]', num2str(chInfo.AveragePathGains, '%.2f '));
    end

    pathPowers = 10.^(chInfo.AveragePathGains/10);
    meanDelay = sum(chInfo.PathDelays .* pathPowers) / sum(pathPowers);
    rmsDelaySpread = sqrt(sum(((chInfo.PathDelays - meanDelay).^2) .* pathPowers) / sum(pathPowers));

    if logging == 1
        logBuffer(end+1) = sprintf('RMS Delay Spread: %.2f ns', rmsDelaySpread * 1e9);
        logBuffer(end+1) = sprintf('Mean Excess Delay: %.2f ns', meanDelay * 1e9);
        if isprop(channel, 'MaximumDopplerShift')
            logBuffer(end+1) = sprintf('Maximum Doppler Shift: %.2f Hz', channel.MaximumDopplerShift);
        end
        logBuffer(end+1) = sprintf('================================================\n');
    end

    %% --------------------------------------------------------------------
    % Transmission / reception setup
    % ---------------------------------------------------------------------
    constPlot = comm.ConstellationDiagram;
    constPlot.EnableMeasurements = 1;

    offset = 0;

    totalInitialBits = 0;
    totalTxBits = 0;
    totalRxBits = 0;
    totalTransmissions = 0;
    totalInitialTransmissions = 0;
    totalRetransmissions = 0;
    totalCodedBits = 0;
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

    % Per-HARQ-process modulation memory:
    processModulation = strings(NHARQProcesses,1);
    processModulation(:) = fixedModulation;

    % Simple feedback feature: recent NACK rate over first attempts only
    feedbackWindow = 50;
    recentFirstAttemptOutcomes = nan(feedbackWindow,1);   % 1=success, 0=need retx
    feedbackPtr = 1;

    tbRows = struct( ...
        'Slot', {}, ...
        'HARQProcessID', {}, ...
        'SNRdB', {}, ...
        'Modulation', {}, ...
        'Is64QAM', {}, ...
        'PredictedRisk', {}, ...
        'RecentNACKRate', {}, ...
        'NeedRetx', {}, ...
        'SuccessfulDecode', {}, ...
        'ThroughputBits', {}, ...
        'TBSize', {}, ...
        'G', {});

    lastPdschEq = [];
    lastRefConst = [];
    selectedRisk = nan;
    recentNACKRate = getRecentNACKRate(recentFirstAttemptOutcomes);

    %% --------------------------------------------------------------------
    % Main simulation loop
    % ---------------------------------------------------------------------
    for nSlot = 0:totalNoSlots-1
        carrier.NSlot = nSlot;
        currentHARQ = harqEntity.HARQProcessID + 1;

        % Decide modulation only when a new TB starts on this HARQ process
        isNewDataCW = false(pdsch.NumCodewords,1);
        for cwIdx = 1:pdsch.NumCodewords
            isNewDataCW(cwIdx) = harqEntity.NewData(cwIdx);
        end

        if any(isNewDataCW)
            recentNACKRate = getRecentNACKRate(recentFirstAttemptOutcomes);

            if controlMode == "ml"
                p16 = predictRetransmissionRisk(mlModel, SNRdB, 0, recentNACKRate);
                p64 = predictRetransmissionRisk(mlModel, SNRdB, 1, recentNACKRate);

                if p64 <= mlRiskThreshold
                    selectedMod = "64QAM";
                elseif p16 <= p64
                    selectedMod = "16QAM";
                else
                    selectedMod = "64QAM";
                end
            else
                selectedMod = fixedModulation;
                p16 = nan; p64 = nan;
            end

            processModulation(currentHARQ) = selectedMod;
        end

        pdsch.Modulation = char(processModulation(currentHARQ));
        constPlot.ReferenceConstellation = getConstellationRefPoints(pdsch.Modulation);

        if strcmpi(pdsch.Modulation,'64QAM')
            selectedRisk = p64;
        else
            selectedRisk = p16;
        end

        [pdschIndices,pdschInfo] = nrPDSCHIndices(carrier,pdsch);
        Xoh_PDSCH = 0;
        trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet),pdschInfo.NREPerPRB,codeRate,Xoh_PDSCH);

        for cwIdx = 1:pdsch.NumCodewords
            if harqEntity.NewData(cwIdx)
                trBlk = randi([0 1],trBlkSizes(cwIdx),1);
                setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);
                totalInitialBits = totalInitialBits + trBlkSizes(cwIdx);

                if harqEntity.SequenceTimeout(cwIdx)
                    resetSoftBuffer(decodeDLSCH,cwIdx-1,harqEntity.HARQProcessID);
                    blockErrors = blockErrors + 1;
                end
            end
        end

        codedTrBlock = encodeDLSCH(pdsch.Modulation,pdsch.NumLayers,pdschInfo.G,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
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
            pdschGrid(dmrsAntIndices) = pdschGrid(dmrsAntIndices) + dmrsSymbols(:,p)*precodingWeights(p,:);
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
            [estChGridLayers,noiseEst] = nrChannelEstimate(carrier,rxGrid,dmrsIndices,dmrsSymbols,'CDMLengths',pdsch.DMRS.CDMLengths);
            estChGridAnts = precodeChannelEstimate(estChGridLayers,conj(precodingWeights));
            newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChGridAnts);
        end

        channelMagnitude = abs(estChGridLayers);
        avgChannelGain(nSlot+1) = 10*log10(mean(channelMagnitude(:).^2));

        for layerIdx = 1:pdsch.NumLayers
            H_layer = squeeze(estChGridLayers(:,:,layerIdx,:));
            H_matrix = squeeze(H_layer(max(1,round(size(H_layer,1)/2)), :, :));
            if ~isempty(H_matrix) && size(H_matrix, 2) >= pdsch.NumLayers
                conditionNumber(nSlot+1) = cond(H_matrix);
            end
        end

        [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxGrid,estChGridLayers);
        [pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);

        refConst_col = getConstellationRefPoints(pdsch.Modulation);
        refConst_col = refConst_col(:);

        for layerIdx = 1:pdsch.NumLayers
            layerSyms = pdschEq(:, layerIdx);
            layerSyms = layerSyms / sqrt(mean(abs(layerSyms).^2));
            M = evmMerMetricsDD(layerSyms, refConst_col);
            EVMrms_pct(nSlot+1, layerIdx) = M.RmsEVM_pct;
            EVMpk_pct (nSlot+1, layerIdx) = M.PeakEVM_pct;
            EVMrms_dB (nSlot+1, layerIdx) = M.AvgEVM_dB;
            EVMpk_dB  (nSlot+1, layerIdx) = M.PeakEVM_dB;
            MERavg_dB (nSlot+1, layerIdx) = M.AvgMER_dB;
        end

        for layerIdx = 1:pdsch.NumLayers
            layerSymbols = pdschEq(:, layerIdx);
            [~, symbolDecisions] = min(abs(layerSymbols - refConst_col.'), [], 2);
            decidedSymbols = refConst_col(symbolDecisions);
            signalPower = mean(abs(decidedSymbols).^2);
            noisePower  = mean(abs(layerSymbols - decidedSymbols).^2);
            instantaneousSINR(nSlot+1, layerIdx) = 10*log10(signalPower / max(noisePower, eps));
        end

        if dispPlots
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
            harqEntity.RedundancyVersion,harqEntity.HARQProcessID);

        % Accounting
        for cwIdx = 1:pdsch.NumCodewords
            totalTransmissions = totalTransmissions + 1;
            totalCodedBits = totalCodedBits + pdschInfo.G;
            totalTxBits = totalTxBits + trBlkSizes(cwIdx);

            if isNewDataCW(cwIdx)
                totalInitialTransmissions = totalInitialTransmissions + 1;
            else
                totalRetransmissions = totalRetransmissions + 1;
            end

            if ~blkerr(cwIdx)
                successfulBlocks = successfulBlocks + 1;
                totalRxBits = totalRxBits + trBlkSizes(cwIdx);
                perSlotSuccess(nSlot+1) = 1;
            end
        end

        % Dataset rows: first-attempt TB only
        for cwIdx = 1:pdsch.NumCodewords
            if isNewDataCW(cwIdx)
                recentNACKRate = getRecentNACKRate(recentFirstAttemptOutcomes);

                newRow.Slot = nSlot;
                newRow.HARQProcessID = harqEntity.HARQProcessID;
                newRow.SNRdB = SNRdB;
                newRow.Modulation = string(pdsch.Modulation);
                newRow.Is64QAM = double(strcmpi(pdsch.Modulation,'64QAM'));
                if strcmpi(pdsch.Modulation,'64QAM')
                    newRow.PredictedRisk = p64;
                else
                    newRow.PredictedRisk = p16;
                end
                newRow.RecentNACKRate = recentNACKRate;
                newRow.NeedRetx = double(blkerr(cwIdx) ~= 0);
                newRow.SuccessfulDecode = double(blkerr(cwIdx) == 0);
                newRow.ThroughputBits = double(~blkerr(cwIdx)) * trBlkSizes(cwIdx);
                newRow.TBSize = trBlkSizes(cwIdx);
                newRow.G = pdschInfo.G;
                tbRows(end+1) = newRow; %#ok<AGROW>

                recentFirstAttemptOutcomes(feedbackPtr) = double(blkerr(cwIdx) == 0);
                feedbackPtr = feedbackPtr + 1;
                if feedbackPtr > feedbackWindow; feedbackPtr = 1; end
            end
        end

        statusReport = updateAndAdvance(harqEntity,blkerr,trBlkSizes,pdschInfo.G);
        [perfState, perfLog(nSlot+1)] = logPerfMetrics(perfState, nSlot, trBlkSizes, blkerr, harqEntity, slotDuration_s);

        perSlotBER(nSlot+1) = sum(blkerr) / pdsch.NumCodewords;

        if logging == 1
            for layerIdx = 1:pdsch.NumLayers
                logBuffer(end+1) = sprintf(['Slot %3d | Layer %d | RMS EVM = %6.2f%% | Peak EVM = %6.2f%% | ' ...
                     'Avg EVM = %6.2f dB | Peak EVM = %6.2f dB | Avg MER = %6.2f dB'], ...
                     nSlot, layerIdx, EVMrms_pct(nSlot+1, layerIdx), EVMpk_pct(nSlot+1, layerIdx), ...
                     EVMrms_dB(nSlot+1, layerIdx), EVMpk_dB(nSlot+1, layerIdx), MERavg_dB(nSlot+1, layerIdx));
            end

            logBuffer(end+1) = sprintf(['Slot %3d | NewTB=%d | TxBits=%6d | RxBits=%6d | ' ...
                 'AttemptBLER=%.3f | CumAttemptBLER=%.3f | FinalBLER=%.3f | Eff=%.1f%% | InstThr=%.2f | AvgThr=%.2f'], ...
                nSlot, perfLog(nSlot+1).isNewTB, perfLog(nSlot+1).txBitsThisSlot, ...
                perfLog(nSlot+1).rxBitsThisSlot, ...
                perfLog(nSlot+1).attemptBLER_thisSlot, perfLog(nSlot+1).attemptBLER_cum, ...
                perfLog(nSlot+1).finalBLER, perfLog(nSlot+1).throughputEfficiencyPct, ...
                perfLog(nSlot+1).instThroughputMbps, perfLog(nSlot+1).avgThroughputMbps);

            if controlMode == "ml"
                logBuffer(end+1) = sprintf('Slot %3d | ML Control | HARQ=%d | Mod=%s | PredRisk=%.4f | RecentNACK=%.4f', ...
                    nSlot, currentHARQ-1, string(pdsch.Modulation), selectedRisk, recentNACKRate);
            end

            logBuffer(end+1) = sprintf('Slot %d. %s', nSlot, statusReport);
        end

        lastPdschEq = pdschEq;
        lastRefConst = refConst_col;
    end

    %% --------------------------------------------------------------------
    % Save ML dataset for this run
    % ---------------------------------------------------------------------
    if datasetMode && ~isempty(tbRows) && strlength(result.datasetFile) > 0
        datasetTable = struct2table(tbRows);
        save(datasetFile, 'datasetTable');
        if logging == 1
            logBuffer(end+1) = sprintf('Dataset file saved to: %s', datasetFile);
        end
    end

    %% --------------------------------------------------------------------
    % Save constellation figure
    % ---------------------------------------------------------------------
    if ~isempty(lastPdschEq)
        refConst = lastRefConst(:);
        layerColors = [0.15 0.45 0.85;
                       0.85 0.33 0.10;
                       0.47 0.67 0.19;
                       0.64 0.08 0.18];

        figConst = figure('Name','Constellation','Visible','off');
        figConst.Position = [100 100 580 560];
        set(figConst, 'Color', 'white');
        ax = axes('Parent', figConst);
        hold(ax, 'on');

        for lyr = 1:pdsch.NumLayers
            layerSyms = lastPdschEq(:, lyr);
            c = layerColors(lyr, :);
            plot(ax, real(layerSyms), imag(layerSyms), '.', ...
                 'Color', c, 'MarkerSize', 3, ...
                 'DisplayName', sprintf('Layer %d Rx', lyr));
        end

        plot(ax, real(refConst), imag(refConst), 'r+', ...
             'MarkerSize', 5, 'LineWidth', 0.8, 'DisplayName', 'Reference');

        hold(ax, 'off');
        axis(ax, 'equal');
        grid(ax, 'on');
        xlabel(ax, 'In-phase');
        ylabel(ax, 'Quadrature');
        title(ax, sprintf('Equalized PDSCH Constellation | Last Slot | SNR = %.1f dB', SNRdB));
        legend(ax, 'Location', 'northeast', 'FontSize', 8);

        if logging == 1
            exportgraphics(figConst, figFileConst, 'Resolution', 150);
            logBuffer(end+1) = sprintf('Constellation figure saved to: %s', figFileConst);
        end
        close(figConst);
    end

    %% --------------------------------------------------------------------
    % Summary logging
    % ---------------------------------------------------------------------
    attemptBLER = perfState.attemptsFailed / max(1, perfState.attemptsTotal);
    finalBLER = 1 - (successfulBlocks / max(1,totalInitialTransmissions));
    throughputEfficiency = totalRxBits / max(1,totalInitialBits);
    harqEfficiency = totalRxBits / max(1,totalTxBits);
    avgBitsPerSlot = totalRxBits / totalNoSlots;
    throughput_Mbps = (avgBitsPerSlot / slotDuration_s) / 1e6;
    latency_est = (1 + (totalRetransmissions / max(1,totalInitialTransmissions))) * slotDuration_s;

    if isempty(tbRows)
        mod16Count = 0;
        mod64Count = 0;
        avgPredRisk = nan;
        avgRecentNack = nan;
        nominalCodeRate = nan;
    else
        tbMods = string({tbRows.Modulation});
        mod16Count = sum(tbMods == "16QAM");
        mod64Count = sum(tbMods == "64QAM");
        avgPredRisk = mean([tbRows.PredictedRisk]);
        avgRecentNack = mean([tbRows.RecentNACKRate]);
        nominalCodeRate = tbRows(1).TBSize / max(1, tbRows(1).G);
    end

    logBuffer(end+1) = sprintf('========================================');
    logBuffer(end+1) = sprintf('     5G NR HARQ Simulation Results      ');
    logBuffer(end+1) = sprintf('========================================');

    logBuffer(end+1) = sprintf('--- Configuration ---');
    if controlMode == "fixed"
        summaryMod = string(fixedModulation);
    else
        summaryMod = "ML_16QAM_64QAM";
    end
    logBuffer(end+1) = sprintf('HARQ Type: Chase Combining');
    logBuffer(end+1) = sprintf('RV Sequence: [%s]', strjoin(string(rvSeq), ' '));
    logBuffer(end+1) = sprintf('Number of HARQ Processes: %d', NHARQProcesses);
    logBuffer(end+1) = sprintf('Modulation: %s', summaryMod);
    logBuffer(end+1) = sprintf('Number of Layers: %d', pdsch.NumLayers);
    logBuffer(end+1) = sprintf('Code Rate: %.4f', codeRate(1));
    logBuffer(end+1) = sprintf('SNR: %.2f dB', SNRdB);
    logBuffer(end+1) = sprintf('Number of Slots: %d', totalNoSlots);
    logBuffer(end+1) = sprintf('Control Mode: %s', controlMode);
    if controlMode == "fixed"
        logBuffer(end+1) = sprintf('Fixed Modulation: %s', fixedModulation);
    else
        logBuffer(end+1) = sprintf('ML Threshold: %.3f', mlRiskThreshold);
        logBuffer(end+1) = sprintf('16QAM TB Count: %d', mod16Count);
        logBuffer(end+1) = sprintf('64QAM TB Count: %d', mod64Count);
        logBuffer(end+1) = sprintf('Average Predicted Risk: %.4f', avgPredRisk);
        logBuffer(end+1) = sprintf('Average Recent NACK Rate: %.4f', avgRecentNack);
    end
    logBuffer(end+1) = sprintf('Dataset Mode: %d', datasetMode);
    if datasetMode
        logBuffer(end+1) = sprintf('Dataset File: %s', datasetFile);
    end

    logBuffer(end+1) = sprintf('--- Transmission Statistics ---');
    logBuffer(end+1) = sprintf('Total Transmissions: %d', totalTransmissions);
    logBuffer(end+1) = sprintf('Initial Transmissions: %d', totalInitialTransmissions);
    logBuffer(end+1) = sprintf('Retransmissions: %d', totalRetransmissions);
    logBuffer(end+1) = sprintf('Average Transmissions per TB: %.2f', ...
        totalTransmissions / max(1,totalInitialTransmissions));
    logBuffer(end+1) = sprintf('Retransmission Rate: %.2f%%', ...
        100 * (totalRetransmissions / max(1,totalInitialTransmissions)));

    logBuffer(end+1) = sprintf('--- Performance Metrics ---');
    logBuffer(end+1) = sprintf('Successful Blocks: %d / %d', successfulBlocks, totalInitialTransmissions);
    logBuffer(end+1) = sprintf('Failed Blocks (after max retx): %d', blockErrors);
    logBuffer(end+1) = sprintf('Attempt BLER (pre-HARQ): %.4f (%.2f%%)', attemptBLER, attemptBLER*100);
    logBuffer(end+1) = sprintf('Final BLER (post-HARQ):  %.4f (%.2f%%)', finalBLER, finalBLER*100);

    logBuffer(end+1) = sprintf('Total Information Bits Sent (all HARQ attempts): %d', totalTxBits);
    logBuffer(end+1) = sprintf('Total Information Bits Delivered Successfully: %d', totalRxBits);
    logBuffer(end+1) = sprintf('Total Coded Bits Transmitted: %d', totalCodedBits);
    logBuffer(end+1) = sprintf('Throughput Efficiency: %.2f%%', throughputEfficiency * 100);
    logBuffer(end+1) = sprintf('HARQ Efficiency: %.2f%%', (totalRxBits / max(1,totalTxBits)) * 100);
    logBuffer(end+1) = sprintf('Effective Code Rate (InfoBits / CodedBits): %.4f', ...
        totalRxBits / max(1,totalCodedBits));
    logBuffer(end+1) = sprintf('Nominal Code Rate (per TB): %.4f', nominalCodeRate);
    logBuffer(end+1) = sprintf('Average Throughput: %.2f Mbps', throughput_Mbps);
    logBuffer(end+1) = sprintf('Average Bits per Slot: %.0f bits', avgBitsPerSlot);

    logBuffer(end+1) = sprintf('--- Latency Calculation ---');
    logBuffer(end+1) = sprintf('Time Slot Duration (ms): %.6f', slotDuration_s * 1e3);
    logBuffer(end+1) = sprintf('Estimated Latency (ms): %.6f', latency_est * 1e3);

    logBuffer(end+1) = sprintf('--- Spectral Efficiency ---');
    numRBs = numel(pdsch.PRBSet);
    bandwidth_Hz = numRBs * 12 * carrier.SubcarrierSpacing * 1e3;
    spectralEfficiency = (avgBitsPerSlot / slotDuration_s) / bandwidth_Hz;
    logBuffer(end+1) = sprintf('Spectral Efficiency: %.4f bits/s/Hz', spectralEfficiency);

    logBuffer(end+1) = sprintf('========================================');
    logBuffer(end+1) = sprintf('--- Channel Model Details ---');
    logBuffer(end+1) = sprintf('Delay Profile: %s', channel.DelayProfile);
    logBuffer(end+1) = sprintf('RMS Delay Spread: %.2f ns', rmsDelaySpread * 1e9);
    logBuffer(end+1) = sprintf('Mean Excess Delay: %.2f ns', meanDelay * 1e9);
    logBuffer(end+1) = sprintf('Maximum Channel Delay: %d samples (%.2f μs)', ...
        chInfo.MaximumChannelDelay, ...
        chInfo.MaximumChannelDelay / channel.SampleRate * 1e6);

    logBuffer(end+1) = sprintf('--- Channel Quality Statistics ---');
    logBuffer(end+1) = sprintf('Average Channel Gain: %.2f dB', mean(avgChannelGain));
    logBuffer(end+1) = sprintf('Channel Gain Std Dev: %.2f dB', std(avgChannelGain));
    logBuffer(end+1) = sprintf('Average Condition Number: %.2f', mean(conditionNumber(conditionNumber > 0)));

    logBuffer(end+1) = sprintf('--- SINR Statistics (per layer) ---');
    for layerIdx = 1:pdsch.NumLayers
        logBuffer(end+1) = sprintf('Layer %d - Mean SINR: %.2f dB, Std: %.2f dB', ...
            layerIdx, mean(instantaneousSINR(:, layerIdx)), std(instantaneousSINR(:, layerIdx)));
    end

    coherenceBW_Hz = 1 / (5 * rmsDelaySpread);
    logBuffer(end+1) = sprintf('Coherence Bandwidth (50%% corr): %.2f kHz', coherenceBW_Hz / 1e3);

    if isprop(channel, 'MaximumDopplerShift') && channel.MaximumDopplerShift > 0
        coherenceTime_s = 9 / (16 * pi * channel.MaximumDopplerShift);
        logBuffer(end+1) = sprintf('Coherence Time (50%% corr): %.2f ms', coherenceTime_s * 1e3);
    end

    %% --------------------------------------------------------------------
    % Return summary values
    % ---------------------------------------------------------------------
    result.ok = true;
    result.ThroughputEfficiencyPct = 100*throughputEfficiency;
    result.HARQEfficiencyPct = 100*harqEfficiency;
    result.RetransmissionRatePct = 100*(totalRetransmissions / max(1,totalInitialTransmissions));
    result.AttemptBLERPct = 100*attemptBLER;
    result.FinalBLERPct = 100*finalBLER;
    result.AverageThroughputMbps = throughput_Mbps;
    result.EstimatedLatencyMs = latency_est * 1e3;
    if controlMode == "fixed"
        result.Modulation = string(fixedModulation);
    elseif controlMode == "ml"
        result.Modulation = "ML_16QAM_64QAM";
    else
        result.Modulation = "UNKNOWN";
    end

    %% --------------------------------------------------------------------
    % Final log flush
    % ---------------------------------------------------------------------
    if logging == 1 && fid >= 0
        fprintf(fid, '%s\n', logBuffer);
        fprintf(fid, '\nRun End: %s\n', datestr(now));
        elapsedTime = toc(runTimer);
        fprintf(fid, '\nTotal Time: %s\n', formatDuration(elapsedTime));
        fclose(fid);
        fid = -1;
    end

catch ME
    if logging == 1 && fid >= 0
        if ~isempty(logBuffer)
            fprintf(fid, '%s\n', logBuffer);
        end
        fprintf(fid, '\nERROR:\n%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        fclose(fid);
        fid = -1;
    end
    rethrow(ME);
end
end

function p = predictRetransmissionRisk(model, SNRdB, is64QAM, recentNACKRate)
if isempty(model)
    % Conservative fallback if model not supplied
    p = 0.5 + 0.2*double(is64QAM) - 0.03*SNRdB + 0.4*recentNACKRate;
    p = min(max(p,0),1);
    return;
end

x = [SNRdB, is64QAM, recentNACKRate];
z = (x - model.mu) ./ model.sigma;
eta = model.beta(1) + z * model.beta(2:end);
p = 1 ./ (1 + exp(-eta));
end

function r = getRecentNACKRate(outcomes)
valid = ~isnan(outcomes);
if ~any(valid)
    r = 0;
else
    successRate = mean(outcomes(valid));
    r = 1 - successRate;
end
end

function str = formatDuration(seconds)
h = floor(seconds/3600);
m = floor(mod(seconds,3600)/60);
s = mod(seconds,60);
str = sprintf('%02dh:%02dm:%05.2fs', h, m, s);
end
