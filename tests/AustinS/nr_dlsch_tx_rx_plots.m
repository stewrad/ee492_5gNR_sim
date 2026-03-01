% DL-SCH and PDSCH Transmit and Receive Processing Chain: 
% This example shows how to use 5G Toolbox™ features to model a 5G NR physical downlink shared channel (PDSCH) link, including all of the steps from transport block generation to bit decoding at the receiver end.

% SET TEST PARAMETERS HERE: 
% ------------------------------------------------------------------------------------------------
SNRdB = 8;                                                                                        % [0 2 4 6 8 10 12] dB // baseline 8.0dB; 
Modulation = "QPSK";                                                                             % must be 'QPSK', '16QAM', '64QAM', '256QAM', '1024QAM'; // baseline QPSK
NHARQProcesses = 16;                                                                               % [1 2 4 8 16]; // baseline is 16
rvSeq = [0 2 3 1];                                                                                 % [0 2 3 1], [0 1 2 3] (monotonic), [0 3 2 1] (reverse), [0 2 1 3]; // originally set as [0 2 3 1]
% Specify the number of transmit and receive antennas 
nTxAnts = 8;                                                                                       % [1 2 4 8] // baseline is Tx = 8 
nRxAnts = 8;                                                                                       % [1 2 4 8] // baseline is Rx = 8
NumLayers = 2;                                                                                     % Set to 2 for all but Tx/Rx Ant = 1;
% Specify Channel Model
DelayProfile = "TDL-C";                                                                            % TDL-A,B,C and CDL-A,B,C // TDL-C baseline
% DL-SCH and PDSCH Transmit and Receive Processing Chain: 
% This example shows how to use 5G Toolbox™ features to model a 5G NR physical downlink shared channel (PDSCH) link, including all of the steps from transport block generation to bit decoding at the receiver end.

% OPTIONAL PLOTTING VARIABLE:
dispPlots = 0;

% OPTIONAL LOGGING VARIABLE: 
logging = 1;

if logging == 1
    % Define log directory name and location, create if it doesn't exist
    logDir = fullfile(pwd, 'logs');  
    if ~exist(logDir,'dir')
        mkdir(logDir);
    end

    % Define figures directory name and location, create if it doesn't exist
    figDir = fullfile(pwd,'figures');
    if ~exist(figDir,'dir')
        mkdir(figDir);
    end

    % Set rvSeq and number of antennas to strings so they can be read in easily as a variable
    rvStr  = sprintf('%d', rvSeq);
    antStr = sprintf('%dx%d', nTxAnts, nRxAnts);
    
    % logFile base naming convention generation 
    logFile = fullfile(logDir, ...
        sprintf("%s_SNR%.1f_%s_NHARQ%d_rvSeq[%s]_%s_%s.txt", ...
        Modulation, ...
        SNRdB, ...
        DelayProfile, ...
        NHARQProcesses, ...
        rvStr, ...
        antStr, ...
        datestr(now,'yyyymmdd_HHMMSS')));

    % Create figures .png with the same base filename as logFile 
    [~, baseName, ~] = fileparts(logFile);
    figFile        = fullfile(figDir, baseName + ".png");           % kept for compatibility
    figFileChEst   = fullfile(figDir, baseName + "_ChEst.png");
    figFileConst   = fullfile(figDir, baseName + "_Constellation.png");

    diary(logFile);
end

% Define a struct of parameters set in this simulation run 
% ------------------------------------------------------------------------------------------------
runParams.SNRdB = SNRdB;
runParams.Modulation = Modulation;
runParams.NHARQProcesses = NHARQProcesses;
runParams.rvSeq = rvSeq;
runParams.nTxAnts = nTxAnts;
runParams.nRxAnts = nRxAnts;
runParams.NumLayers = NumLayers;
runParams.DelayProfile = DelayProfile;
% runParams.seed = seed;
% [~, gitHash] = system('git rev-parse HEAD');



% Specify SNR, number of slots, and perfect channel estimation flag 
% ------------------------------------------------------------------------------------------------
SNRdB = SNRdB;             % SNR in dB
totalNoSlots = 20;         % Number of slots to simulate
perfectEstimation = false; % Perfect synchronization and channel estimation
rng("default");            % Set default random number generator for repeatability

% Carrier configuration:
% Create config object. this object controls the numerology, i.e. subcarrier spacing, bandwidth, and cyclic prefix (CP) length. 
% ------------------------------------------------------------------------------------------------
carrier = nrCarrierConfig;

% PDSCH and DM-RS config: 
% Create a PDSCH config object. Specify modulation scheme, number of layers.. allocate all resource blocks (RBs) to PDSCH (full band allocation) 
% ------------------------------------------------------------------------------------------------
pdsch = nrPDSCHConfig;
pdsch.Modulation = Modulation;
pdsch.NumLayers = NumLayers;
pdsch.PRBSet = 0:carrier.NSizeGrid-1;     % Full band allocation

% Set DM-RS parameters 
% Add additional DM-RS position to improve channel estimation 
% DMRSConfigurationType = 1 supports up to 4 DM-RS ports when DMRSLength = 1.
% DMRSConfigurationType = 1 supports up to 8 DM-RS ports when DMRSLength = 2.
% DMRSConfigurationType = 2 supports up to 6 DM-RS ports when DMRSLength = 1. This is designed for multi-user MIMO (MU-MIMO).
% DMRSConfigurationType = 2 supports up to 12 DM-RS ports when DMRSLength = 2. This is designed for MU-MIMO.
% ------------------------------------------------------------------------------------------------
pdsch.DMRS.DMRSAdditionalPosition = 1;
pdsch.DMRS.DMRSConfigurationType = 1;
pdsch.DMRS.DMRSLength = 2;
% pdsch.DMRS;                            % Display DM-RS properties

% DL-SCH Configuration 
% Specify code rate, number of HARQ processes, and redundancy version (RV) sequence values. 
% Note: to disable HARQ retransmissions, set rvSeq to a fix value (e.g. 0) 
% ------------------------------------------------------------------------------------------------
% NHARQProcesses = 16;     % Number of parallel HARQ processes
NHARQProcesses = NHARQProcesses;     % Number of parallel HARQ processes
% rvSeq = [0 2 3 1];
rvSeq = rvSeq; 

% Coding rate
% Take into account # of codewords when specifying code rate: 
% 1 codeword for up to 4 layers
% 2 codewords for more than 4 layers
% ------------------------------------------------------------------------------------------------
if pdsch.NumCodewords == 1
    codeRate = 490/1024;
else
    codeRate = [490 490]./1024;
end

% Create DL-SCH encoder and decoder objects 
% To use multiple processes, set MultipleHARQProcesses property to true for both objects
% DL-SCH encoder/decode objects can handle up to 16 HARQ processes 
% ------------------------------------------------------------------------------------------------
% Create DL-SCH encoder object
encodeDLSCH = nrDLSCH;
encodeDLSCH.MultipleHARQProcesses = true;
encodeDLSCH.TargetCodeRate = codeRate;
% Create DLSCH decoder object
decodeDLSCH = nrDLSCHDecoder;
decodeDLSCH.MultipleHARQProcesses = true;
decodeDLSCH.TargetCodeRate = codeRate;
decodeDLSCH.LDPCDecodingAlgorithm = "Normalized min-sum";
decodeDLSCH.MaximumLDPCIterationCount = 6;

% HARQ Management
% Create a HARQ entity object to manage the HARQ processes and DL-SCH encoder and decoder buffers 
harqEntity = HARQEntity(0:NHARQProcesses-1,rvSeq,pdsch.NumCodewords);

% Channel Configuration
% ------------------------------------------------------------------------------------------------
% Specify the number of transmit and receive antennas 
% nTxAnts = 8;
nTxAnts = nTxAnts;
% nRxAnts = 8;
nRxAnts = nRxAnts;

% Check that the number of layers is valid for the number of antennas
if pdsch.NumLayers > min(nTxAnts,nRxAnts)
    error("The number of layers ("+string(pdsch.NumLayers)+") must be smaller than min(nTxAnts,nRxAnts) ("+string(min(nTxAnts,nRxAnts))+")")
end

% % Create a channel object 
% channel = nrTDLChannel;
% channel.DelayProfile = "TDL-C";
% channel.NumTransmitAntennas = nTxAnts;
% channel.NumReceiveAntennas = nRxAnts;
% Create a channel object 
channel = nrTDLChannel;
channel.DelayProfile = DelayProfile;
channel.NumTransmitAntennas = nTxAnts;
channel.NumReceiveAntennas = nRxAnts;

% Set the channel sample rate of the OFDM signal 
% Note: to obtain the sample rate of the OFDM signal, use nrOFDMInfo function 
ofdmInfo = nrOFDMInfo(carrier);
channel.SampleRate = ofdmInfo.SampleRate;

% Set the channel output type so that perfect channel estimation and timing estimation can be calculated at the same time when filtering the signal 
channel.ChannelResponseOutput = 'ofdm-response';

% ========== Extract Channel Model Information ==========
chInfo = info(channel);

% Print simulated parameters for test log
fprintf('\n===== RUN START =====\n');
fprintf('Timestamp: %s\n', datestr(now));
disp(runParams);
% fprintf('Git Commit: %s\n', strtrim(gitHash));
fprintf('=====================\n\n');

% Get additional channel parameters
fprintf('\n========== Channel Model Configuration ==========\n');
fprintf('Channel Type: %s\n', channel.DelayProfile);
fprintf('Sample Rate: %.2f MHz\n', channel.SampleRate / 1e6);
fprintf('Maximum Channel Delay: %d samples\n', chInfo.MaximumChannelDelay);
fprintf('Transmit Antennas: %d\n', channel.NumTransmitAntennas);
fprintf('Receive Antennas: %d\n', channel.NumReceiveAntennas);

% Display path details
fprintf('\n--- Path Characteristics ---\n');
fprintf('Path Delays (ns): [%s]\n', num2str(chInfo.PathDelays * 1e9, '%.2f '));
fprintf('Average Path Gains (dB): [%s]\n', num2str(chInfo.AveragePathGains, '%.2f '));

% Calculate RMS delay spread
pathPowers = 10.^(chInfo.AveragePathGains/10);
meanDelay = sum(chInfo.PathDelays .* pathPowers) / sum(pathPowers);
rmsDelaySpread = sqrt(sum(((chInfo.PathDelays - meanDelay).^2) .* pathPowers) / sum(pathPowers));

fprintf('RMS Delay Spread: %.2f ns\n', rmsDelaySpread * 1e9);
fprintf('Mean Excess Delay: %.2f ns\n', meanDelay * 1e9);

% Check if Doppler is configured
if isprop(channel, 'MaximumDopplerShift')
    fprintf('Maximum Doppler Shift: %.2f Hz\n', channel.MaximumDopplerShift);
end

fprintf('================================================\n\n');

% Transmission and reception
% Set up a loop to simulate the transmission and reception of slots. Create a comm.ConstellationDiagram to display the constellation of the equalized signal.
% ------------------------------------------------------------------------------------------------
constPlot = comm.ConstellationDiagram;                                          % Constellation diagram object
constPlot.ReferenceConstellation = getConstellationRefPoints(pdsch.Modulation); % Reference constellation values
constPlot.EnableMeasurements = 1;                                               % Enable EVM measurements

% Initial timing offset
offset = 0;

% ========== Performance Tracking Variables ==========
totalTxBits = 0;           % Total bits transmitted
totalRxBits = 0;           % Total bits successfully received
totalTransmissions = 0;    % Total number of transmissions (including retransmissions)
totalInitialTransmissions = 0;  % Count of new transport blocks
totalRetransmissions = 0;  % Count of retransmissions
blockErrors = 0;           % Number of blocks that failed after all retransmissions
successfulBlocks = 0;      % Number of blocks decoded successfully
transmissionAttempts = zeros(totalNoSlots, 1);  % Track attempts per TB

% For tracking per-transmission statistics
perSlotSuccess = zeros(totalNoSlots, 1);
perSlotBER = zeros(totalNoSlots, 1);

estChannelGrid = getInitialChannelEstimate(channel,carrier);
newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChannelGrid);

% ========== Channel Statistics Tracking ==========
instantaneousSINR = zeros(totalNoSlots, pdsch.NumLayers);
avgChannelGain = zeros(totalNoSlots, 1);
conditionNumber = zeros(totalNoSlots, 1);


% ========== Per-Transmission (Per-Slot) Logging Storage ==========
perfLog = repmat(struct( ...
    'slot', 0, ...
    'isNewTB', false, ...
    'txBitsThisSlot', 0, ...
    'rxBitsThisSlot', 0, ...
    'blkErrCountThisSlot', 0, ...
    'totalTransmissions', 0, ...
    'totalInitialTransmissions', 0, ...
    'totalRetransmissions', 0, ...
    'successfulBlocks', 0, ...
    'blockErrors', 0, ...
    'bler', 0, ...
    'throughputEfficiencyPct', 0, ...
    'instThroughputMbps', 0, ...
    'avgThroughputMbps', 0), totalNoSlots, 1);

% Slot duration for throughput calculations (seconds)
% (Keep your carrier config as-is; this uses the numerology already set)
slotDuration_s = 1e-3 / carrier.SlotsPerSubframe;  % 1 subframe = 1 ms


% ========== Cumulative performance state (for per-slot logging) ==========
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

% ========== Per-slot log storage ==========
% ========== Per-Transmission (Per-Slot) Logging Storage ==========
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


% ================= MAIN SIMULATION LOOP =================
for nSlot = 0:totalNoSlots-1
    % New slot
    carrier.NSlot = nSlot;

    % Calculate transport block size 
    % Generate PDSCH indices info, which is needed to calculate the transport
    % block size
    [pdschIndices,pdschInfo] = nrPDSCHIndices(carrier,pdsch);

    % Calculate transport block sizes
    Xoh_PDSCH = 0;
    trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet),pdschInfo.NREPerPRB,codeRate,Xoh_PDSCH);

    % Get new transport blocks and flush decoder soft buffer, as required
    for cwIdx = 1:pdsch.NumCodewords
        if harqEntity.NewData(cwIdx)
            % Create and store a new transport block for transmission
            trBlk = randi([0 1],trBlkSizes(cwIdx),1);
            setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);

            % If the previous RV sequence ends without successful
            % decoding, flush the soft buffer
            if harqEntity.SequenceTimeout(cwIdx)
                resetSoftBuffer(decodeDLSCH,cwIdx-1,harqEntity.HARQProcessID);
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

    % PDSCH DM-RS precoding and mapping
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
        % Perfect timing estimation is provided by the channel
        offset = timingOffset;
    else
        [t,mag] = nrTimingEstimate(carrier,rxWaveform,dmrsIndices,dmrsSymbols);
        offset = hSkipWeakTimingOffset(offset,t,mag);
    end
    rxWaveform = rxWaveform(1+offset:end,:);

    rxGrid = nrOFDMDemodulate(carrier,rxWaveform);

    if perfectEstimation
        % Perfect channel estimation between transmit and receive antennas
        % provided by the channel
        estChGridAnts = ofdmChannelResponse;

        % Use the precalculated noise variance as the perfect noise
        % estimate
        noiseEst = nVar;

        % Get precoding matrix for next slot
        newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChGridAnts);

        % Apply precoding to estChGridAnts. The resulting estimate is for
        % the channel estimate between layers and receive antennas.
        estChGridLayers = precodeChannelEstimate(estChGridAnts,precodingWeights.');
    else
        % Perform practical channel estimation between layers and receive
        % antennas.
        [estChGridLayers,noiseEst] = nrChannelEstimate(carrier,rxGrid,dmrsIndices,dmrsSymbols,'CDMLengths',pdsch.DMRS.CDMLengths);

        % Remove precoding from estChannelGrid before precoding
        % matrix calculation
        estChGridAnts = precodeChannelEstimate(estChGridLayers,conj(precodingWeights));

        % Get precoding matrix for next slot
        newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChGridAnts);
    end

    % ========== Track Channel Characteristics per Slot ==========
    % Calculate average channel gain (in dB)
    channelMagnitude = abs(estChGridLayers);
    avgChannelGain(nSlot+1) = 10*log10(mean(channelMagnitude(:).^2));
    
    % Calculate condition number for MIMO channel health
    for layerIdx = 1:pdsch.NumLayers
        H_layer = squeeze(estChGridLayers(:,:,layerIdx,:));
        % Get a representative subcarrier
        H_matrix = squeeze(H_layer(size(H_layer,1)/2, :, :)); % Middle subcarrier
        if ~isempty(H_matrix) && size(H_matrix, 2) >= pdsch.NumLayers
            conditionNumber(nSlot+1) = cond(H_matrix);
        end
    end
    
    [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxGrid,estChGridLayers);
    [pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);

    % ---- EVM / MER per layer via evmMerMetricsDD ----
    refConst_col = getConstellationRefPoints(pdsch.Modulation);
    refConst_col = refConst_col(:);   % ensure Kx1 column vector

    for layerIdx = 1:pdsch.NumLayers
        layerSyms = pdschEq(:, layerIdx);

        M = evmMerMetricsDD(layerSyms, refConst_col);

        EVMrms_pct(nSlot+1, layerIdx) = M.RmsEVM_pct;
        EVMpk_pct (nSlot+1, layerIdx) = M.PeakEVM_pct;
        EVMrms_dB (nSlot+1, layerIdx) = M.AvgEVM_dB;
        EVMpk_dB  (nSlot+1, layerIdx) = M.PeakEVM_dB;
        MERavg_dB (nSlot+1, layerIdx) = M.AvgMER_dB;

        fprintf(['Slot %3d | Layer %d | RMS EVM = %6.2f%% | Peak EVM = %6.2f%% | ' ...
                 'Avg EVM = %6.2f dB | Peak EVM = %6.2f dB | Avg MER = %6.2f dB\n'], ...
                 nSlot, layerIdx, M.RmsEVM_pct, M.PeakEVM_pct, M.AvgEVM_dB, M.PeakEVM_dB, M.AvgMER_dB);
    end

    % Estimate instantaneous SINR per layer from equalized symbols
    for layerIdx = 1:pdsch.NumLayers
        layerSymbols = pdschEq(:, layerIdx);

        % Find nearest constellation point (reuse refConst_col from EVM block)
        [~, symbolDecisions] = min(abs(layerSymbols - refConst_col.'), [], 2);
        decidedSymbols = refConst_col(symbolDecisions);

        signalPower = mean(abs(decidedSymbols).^2);
        noisePower  = mean(abs(layerSymbols - decidedSymbols).^2);

        instantaneousSINR(nSlot+1, layerIdx) = 10*log10(signalPower / noisePower);
    end

    if dispPlots 
        mesh(abs(estChGridLayers(:,:,1,1)));
        title('Channel Estimate');
        xlabel('OFDM Symbol');
        ylabel("Subcarrier");
        zlabel("Magnitude");
    
        constPlot.ChannelNames = "Layer "+(pdsch.NumLayers:-1:1);
        constPlot.ShowLegend = true;
        % Constellation for the first layer has a higher SNR than that for the
        % last layer. Flip the layers so that the constellations do not mask
        % each other.
        constPlot(fliplr(pdschEq));
    end
    
    [dlschLLRs,rxSymbols] = nrPDSCHDecode(carrier,pdsch,pdschEq,noiseEst);

    % Scale LLRs by CSI
    csi = nrLayerDemap(csi);                                    % CSI layer demapping
    for cwIdx = 1:pdsch.NumCodewords
        Qm = length(dlschLLRs{cwIdx})/length(rxSymbols{cwIdx}); % Bits per symbol
        csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);                 % Expand by each bit per symbol
        dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:);   % Scale
    end

    decodeDLSCH.TransportBlockLength = trBlkSizes;
    [decbits,blkerr] = decodeDLSCH(dlschLLRs,pdsch.Modulation,pdsch.NumLayers, ...
        harqEntity.RedundancyVersion,harqEntity.HARQProcessID);

    statusReport = updateAndAdvance(harqEntity,blkerr,trBlkSizes,pdschInfo.G);    

    % ===== Per-transmission (per-slot) metrics snapshot =====
    [perfState, perfLog(nSlot+1)] = logPerfMetrics(perfState, nSlot, trBlkSizes, blkerr, harqEntity, slotDuration_s);

    % Print a compact line every slot (adjust formatting however you want)
    fprintf(['Slot %3d | NewTB=%d | TxBits=%6d | RxBits=%6d | ' ...
         'AttemptBLER=%.3f | CumAttemptBLER=%.3f | ' ...
         'FinalBLER=%.3f | Eff=%.1f%% | InstThr=%.2f | AvgThr=%.2f\n'], ...
        nSlot, perfLog(nSlot+1).isNewTB, perfLog(nSlot+1).txBitsThisSlot, ...
        perfLog(nSlot+1).rxBitsThisSlot, ...
        perfLog(nSlot+1).attemptBLER_thisSlot, perfLog(nSlot+1).attemptBLER_cum, ...
        perfLog(nSlot+1).finalBLER, perfLog(nSlot+1).throughputEfficiencyPct, ...
        perfLog(nSlot+1).instThroughputMbps, perfLog(nSlot+1).avgThroughputMbps);


    % ========== Track Performance Metrics ==========
    totalTransmissions = totalTransmissions + 1;

    for cwIdx = 1:pdsch.NumCodewords
        if harqEntity.NewData(cwIdx)
            totalInitialTransmissions = totalInitialTransmissions + 1;
            totalTxBits = totalTxBits + trBlkSizes(cwIdx);
        else
            totalRetransmissions = totalRetransmissions + 1;
        end

        % Track successful decoding
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

    % Calculate instantaneous BER (if you store original TB)
    % Note: This requires storing the original transport block for comparison
    % For now, we'll track based on block errors
    perSlotBER(nSlot+1) = sum(blkerr) / pdsch.NumCodewords;

    disp("Slot "+(nSlot)+". "+statusReport);

    % Cache last-slot data for end-of-run plotting
    lastEstChGridLayers = estChGridLayers;
    lastPdschEq         = pdschEq;
end % for nSlot = 0:totalNoSlots




% ========== Save Channel Estimate and Constellation Figures ==========
% --- 1. Channel Estimate (mesh of last slot) ---
figChEst = figure('Name','Channel Estimate','Visible','off');
mesh(abs(lastEstChGridLayers(:,:,1,1)));
title('Channel Estimate - Last Slot');
xlabel('OFDM Symbol');
ylabel('Subcarrier');
zlabel('Magnitude');
colorbar;

if logging == 1
    exportgraphics(figChEst, figFileChEst, 'Resolution', 150);
    fprintf('Channel estimate figure saved to: %s\n', figFileChEst);
else
    exportgraphics(figChEst, fullfile(pwd, 'Channel_Estimate.png'), 'Resolution', 150);
    fprintf('Channel estimate figure saved to: %s\n', fullfile(pwd,'Channel_Estimate.png'));
end
close(figChEst);

% --- 2. All-layer Constellation on a single white-background axes ---
refConst   = getConstellationRefPoints(pdsch.Modulation);
refConst   = refConst(:);
layerColors = [0.15 0.45 0.85;   % Layer 1 – blue
               0.85 0.33 0.10;   % Layer 2 – orange
               0.47 0.67 0.19;   % Layer 3 – green  (if ever used)
               0.64 0.08 0.18];  % Layer 4 – dark red

figConst = figure('Name','Constellation','Visible','off');
figConst.Position = [100 100 580 560];

% White figure and axes backgrounds
set(figConst, 'Color', 'white');
ax = axes('Parent', figConst);
set(ax, 'Color', 'white', 'XColor', 'black', 'YColor', 'black', ...
        'GridColor', [0.82 0.82 0.82], 'FontSize', 9);
hold(ax, 'on');

% Overlay each layer with a distinct colour (plotted first so ref sits on top)
for lyr = 1:pdsch.NumLayers
    layerSyms = lastPdschEq(:, lyr);
    c = layerColors(lyr, :);
    plot(ax, real(layerSyms), imag(layerSyms), '.', ...
         'Color', c, 'MarkerSize', 3, ...
         'DisplayName', sprintf('Layer %d Rx', lyr));
end

% Reference symbol positions as red '+' markers (on top)
plot(ax, real(refConst), imag(refConst), 'r+', ...
     'MarkerSize', 5, 'LineWidth', 0.8, 'DisplayName', 'Reference');

hold(ax, 'off');
axis(ax, 'equal');
grid(ax, 'on');
xlabel(ax, 'In-phase Amplitude',    'Color', 'black');
ylabel(ax, 'Quadrature Amplitude',  'Color', 'black');
title(ax, sprintf('Equalized PDSCH Constellation  |  %s  |  SNR = %.1f dB  |  Last Slot', ...
      pdsch.Modulation, SNRdB), 'Color', 'black');

% Legend (white background, black text)
lgd = legend(ax, 'Location', 'northeast', 'FontSize', 8);
set(lgd, 'Color', 'white', 'EdgeColor', 'black', 'TextColor', 'black');

% ---- Combined EVM / MER annotation box (all layers, black text on white) ----
% Row order matches comm.ConstellationDiagram Measurements panel:
%   RMS EVM (%)  |  Peak EVM (%)  |  Avg EVM (dB)  |  Peak EVM (dB)  |  Avg MER (dB)
ann_lines = {};
for lyr = 1:pdsch.NumLayers
    rmsEVM  = EVMrms_pct(end, lyr);
    pkEVM   = EVMpk_pct (end, lyr);
    evmDB   = EVMrms_dB (end, lyr);
    pkEvmDB = EVMpk_dB  (end, lyr);
    merAvg  = MERavg_dB (end, lyr);
    ann_lines{end+1} = sprintf('--- Layer %d ---',         lyr);           %#ok<SAGROW>
    ann_lines{end+1} = sprintf('RMS EVM  : %6.2f %%',     rmsEVM);        %#ok<SAGROW>
    ann_lines{end+1} = sprintf('Peak EVM : %6.2f %%',     pkEVM);         %#ok<SAGROW>
    ann_lines{end+1} = sprintf('Avg EVM  : %6.2f dB',     evmDB);         %#ok<SAGROW>
    ann_lines{end+1} = sprintf('Peak EVM : %6.2f dB',     pkEvmDB);       %#ok<SAGROW>
    ann_lines{end+1} = sprintf('Avg MER  : %6.2f dB',     merAvg);        %#ok<SAGROW>
    if lyr < pdsch.NumLayers
        ann_lines{end+1} = '';   % blank separator between layers           %#ok<SAGROW>
    end
end
ann_str = strjoin(ann_lines, '\n');

xlims = xlim(ax); ylims = ylim(ax);
text(ax, xlims(1) + 0.03*diff(xlims), ylims(1) + 0.04*diff(ylims), ann_str, ...
     'FontSize', 8, 'FontName', 'Courier New', ...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
     'Color', 'black', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
     'Margin', 4);

if logging == 1
    exportgraphics(figConst, figFileConst, 'Resolution', 150);
    fprintf('Constellation figure saved to: %s\n', figFileConst);
else
    exportgraphics(figConst, fullfile(pwd, 'Constellation.png'), 'Resolution', 150);
    fprintf('Constellation figure saved to: %s\n', fullfile(pwd,'Constellation.png'));
end
close(figConst);


% ========== Enhanced Results Display ==========
fprintf('\n========================================\n');
fprintf('     5G NR HARQ Simulation Results      \n');
fprintf('========================================\n');

fprintf('\n--- Configuration ---\n');
fprintf('HARQ Type: Chase Combining\n');
fprintf('RV Sequence: [%s]\n', num2str(rvSeq));
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
fprintf('Average Transmissions per TB: %.2f\n', ...
    totalTransmissions / totalInitialTransmissions);
fprintf('Retransmission Rate: %.2f%%\n', ...
    (totalRetransmissions / totalTransmissions) * 100);

% fprintf('\n--- Performance Metrics ---\n');
% fprintf('Successful Blocks: %d / %d\n', successfulBlocks, totalInitialTransmissions);
% fprintf('Failed Blocks (after max retx): %d\n', blockErrors);
% fprintf('Block Error Rate (BLER): %.4f (%.2f%%)\n', ...
%     blockErrors / totalInitialTransmissions, ...
%     (blockErrors / totalInitialTransmissions) * 100);



% --- BLER metrics ---
finalBLER = blockErrors / max(1,totalInitialTransmissions);              % post-HARQ (what you have now)
attemptBLER = perfState.attemptsFailed / max(1, perfState.attemptsTotal); % pre-HARQ (counts all failed attempts)

fprintf('\n--- Performance Metrics ---\n');
fprintf('Successful Blocks: %d / %d\n', successfulBlocks, totalInitialTransmissions);
fprintf('Failed Blocks (after max retx): %d\n', blockErrors);

fprintf('Attempt BLER (pre-HARQ): %.4f (%.2f%%)\n', attemptBLER, attemptBLER*100);
fprintf('Final BLER (post-HARQ):  %.4f (%.2f%%)\n', finalBLER, finalBLER*100);



fprintf('\n--- Throughput Analysis ---\n');
fprintf('Total Bits Transmitted: %d\n', totalTxBits);
fprintf('Total Bits Received (success): %d\n', totalRxBits);
fprintf('Throughput Efficiency: %.2f%%\n', (totalRxBits / totalTxBits) * 100);
fprintf('Effective Code Rate: %.4f\n', totalRxBits / (totalTransmissions * pdschInfo.G));

% Calculate average throughput per slot
avgBitsPerSlot = totalRxBits / totalNoSlots;
slotDuration = carrier.SlotsPerSubframe / (carrier.SubcarrierSpacing / 15e3) * 1e-3; % seconds
throughput_Mbps = (avgBitsPerSlot / slotDuration) / 1e6;

fprintf('Average Throughput: %.2f Mbps\n', throughput_Mbps);
fprintf('Average Bits per Slot: %.0f bits\n', avgBitsPerSlot);

fprintf('\n--- Spectral Efficiency ---\n');
numRBs = numel(pdsch.PRBSet);
bandwidth_Hz = numRBs * 12 * carrier.SubcarrierSpacing * 1e3;
spectralEfficiency = (avgBitsPerSlot / slotDuration) / bandwidth_Hz;
fprintf('Spectral Efficiency: %.4f bits/s/Hz\n', spectralEfficiency);

fprintf('\n========================================\n');

fprintf('\n--- Channel Model Details ---\n');
fprintf('Delay Profile: %s\n', channel.DelayProfile);
fprintf('RMS Delay Spread: %.2f ns\n', rmsDelaySpread * 1e9);
fprintf('Mean Excess Delay: %.2f ns\n', meanDelay * 1e9);
fprintf('Maximum Channel Delay: %d samples (%.2f μs)\n', ...
    chInfo.MaximumChannelDelay, ...
    chInfo.MaximumChannelDelay / channel.SampleRate * 1e6);

fprintf('\n--- Channel Quality Statistics ---\n');
fprintf('Average Channel Gain: %.2f dB\n', mean(avgChannelGain));
fprintf('Channel Gain Std Dev: %.2f dB\n', std(avgChannelGain));
fprintf('Average Condition Number: %.2f\n', mean(conditionNumber(conditionNumber > 0)));

fprintf('\n--- SINR Statistics (per layer) ---\n');
for layerIdx = 1:pdsch.NumLayers
    fprintf('Layer %d - Mean SINR: %.2f dB, Std: %.2f dB\n', ...
        layerIdx, ...
        mean(instantaneousSINR(:, layerIdx)), ...
        std(instantaneousSINR(:, layerIdx)));
end

% Calculate coherence bandwidth and time
coherenceBW_Hz = 1 / (5 * rmsDelaySpread); % 50% coherence
fprintf('\nCoherence Bandwidth (50%% corr): %.2f kHz\n', coherenceBW_Hz / 1e3);

if isprop(channel, 'MaximumDopplerShift') && channel.MaximumDopplerShift > 0
    coherenceTime_s = 9 / (16 * pi * channel.MaximumDopplerShift);
    fprintf('Coherence Time (50%% corr): %.2f ms\n', coherenceTime_s * 1e3);
end


% Final line for Command Window tracking
if logging == 1
    diary off;
end

