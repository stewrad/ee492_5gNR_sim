% DL-SCH and PDSCH Transmit and Receive Processing Chain: 
% This example shows how to use 5G Toolbox™ features to model a 5G NR physical downlink shared channel (PDSCH) link, including all of the steps from transport block generation to bit decoding at the receiver end.

% SET TEST PARAMETERS HERE: 
% ------------------------------------------------------------------------------------------------
SNRdB = 0;                                                                                        % baseline is 8.0 dB
Modulation = "256QAM";                                                                            % must be 'QPSK', '16QAM', '64QAM', '256QAM', '1024QAM'
NHARQProcesses = 0;                                                                              % baseline is 2
rvSeq = [0 2 3 1];                                                                                % originally set as [0 2 3 1]
nTxAnts = 8;                                                                                      % baseline is Tx/Rx = 8 
nRxAnts = 8;           
NumLayers = 2;                                                                                    % Set to 2 for all but Tx/Rx Ant = 1
DelayProfile = "TDL-C";

% OPTIONAL LOGGING VARIABLE: 
logging = 1;

% HARQ toggle
disableHARQ = true;   % true = NO HARQ, false = normal HARQ

if logging == 1
    logDir = fullfile(pwd, 'logs');  % <-- use lowercase logs
    if ~exist(logDir,'dir')
        mkdir(logDir);
    end
    logFile = fullfile(logDir, sprintf('%s_SNR%.1f_%s_%s.txt', Modulation, SNRdB, DelayProfile, datestr(now,'yyyymmdd_HHMMSS')));
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
runParams.disableHARQ = disableHARQ;
% runParams.seed = seed;
% [~, gitHash] = system('git rev-parse HEAD');

% Specify SNR, number of slots, and perfect channel estimation flag 
% ------------------------------------------------------------------------------------------------
SNRdB = SNRdB;             % SNR in dB
totalNoSlots = 50;         % Number of slots to simulate
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

% DL-SCH Configuration 
% Specify code rate, number of HARQ processes, and redundancy version (RV) sequence values. 
% Note: to disable HARQ retransmissions, set rvSeq to a fixed value (e.g. 0) 
% ------------------------------------------------------------------------------------------------
if disableHARQ
    NHARQProcesses = 1;   % must be at least 1 so HARQEntity can still be created
    rvSeq = 0;            % fixed RV = no retransmissions
end

NHARQProcesses = NHARQProcesses;
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
encodeDLSCH = nrDLSCH;
encodeDLSCH.MultipleHARQProcesses = true;
encodeDLSCH.TargetCodeRate = codeRate;

decodeDLSCH = nrDLSCHDecoder;
decodeDLSCH.MultipleHARQProcesses = true;
decodeDLSCH.TargetCodeRate = codeRate;
decodeDLSCH.LDPCDecodingAlgorithm = "Normalized min-sum";
decodeDLSCH.MaximumLDPCIterationCount = 6;

% HARQ Management
harqEntity = HARQEntity(0:NHARQProcesses-1,rvSeq,pdsch.NumCodewords);

% Channel Configuration
% ------------------------------------------------------------------------------------------------
nTxAnts = nTxAnts;
nRxAnts = nRxAnts;

if pdsch.NumLayers > min(nTxAnts,nRxAnts)
    error("The number of layers ("+string(pdsch.NumLayers)+") must be smaller than min(nTxAnts,nRxAnts) ("+string(min(nTxAnts,nRxAnts))+")")
end

channel = nrTDLChannel;
channel.DelayProfile = DelayProfile;
channel.NumTransmitAntennas = nTxAnts;
channel.NumReceiveAntennas = nRxAnts;

ofdmInfo = nrOFDMInfo(carrier);
channel.SampleRate = ofdmInfo.SampleRate;
channel.ChannelResponseOutput = 'ofdm-response';

% ========== Extract Channel Model Information ==========
chInfo = info(channel);

% Print simulated parameters for test log   
fprintf('\n===== RUN START =====\n');
fprintf('Timestamp: %s\n', datestr(now));
disp(runParams);
fprintf('=====================\n\n');

% Get additional channel parameters
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

% Transmission and reception
% ------------------------------------------------------------------------------------------------
constPlot = comm.ConstellationDiagram;
constPlot.ReferenceConstellation = getConstellationRefPoints(pdsch.Modulation);
constPlot.EnableMeasurements = 1;

offset = 0;

% ========== Performance Tracking Variables ==========
totalTxBits = 0;
totalRxBits = 0;
totalTransmissions = 0;
totalInitialTransmissions = 0;
totalRetransmissions = 0;
blockErrors = 0;
successfulBlocks = 0;
transmissionAttempts = zeros(totalNoSlots, 1);

perSlotSuccess = zeros(totalNoSlots, 1);
perSlotBER = zeros(totalNoSlots, 1);

estChannelGrid = getInitialChannelEstimate(channel,carrier);
newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChannelGrid);

% ========== Channel Statistics Tracking ==========
instantaneousSINR = zeros(totalNoSlots, pdsch.NumLayers);
avgChannelGain = zeros(totalNoSlots, 1);
conditionNumber = zeros(totalNoSlots, 1);

for nSlot = 0:totalNoSlots-1
    carrier.NSlot = nSlot;

    [pdschIndices,pdschInfo] = nrPDSCHIndices(carrier,pdsch);

    Xoh_PDSCH = 0;
    trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet),pdschInfo.NREPerPRB,codeRate,Xoh_PDSCH);

    % Get new transport blocks and flush decoder soft buffer, as required
    for cwIdx = 1:pdsch.NumCodewords
        if disableHARQ || harqEntity.NewData(cwIdx)
            trBlk = randi([0 1],trBlkSizes(cwIdx),1);
            setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);

            if ~disableHARQ && harqEntity.SequenceTimeout(cwIdx)
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

    % ========== Track Channel Characteristics per Slot ==========
    channelMagnitude = abs(estChGridLayers);
    avgChannelGain(nSlot+1) = 10*log10(mean(channelMagnitude(:).^2));
    
    for layerIdx = 1:pdsch.NumLayers
        H_layer = squeeze(estChGridLayers(:,:,layerIdx,:));
        H_matrix = squeeze(H_layer(size(H_layer,1)/2, :, :));
        if ~isempty(H_matrix) && size(H_matrix, 2) >= pdsch.NumLayers
            conditionNumber(nSlot+1) = cond(H_matrix);
        end
    end
    
    [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxGrid,estChGridLayers);
    [pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);

    for layerIdx = 1:pdsch.NumLayers
        layerSymbols = pdschEq(:, layerIdx);
        refConstellation = getConstellationRefPoints(pdsch.Modulation);
        
        [~, symbolDecisions] = min(abs(layerSymbols - refConstellation.'), [], 2);
        decidedSymbols = refConstellation(symbolDecisions);
        
        signalPower = mean(abs(decidedSymbols).^2);
        noisePower = mean(abs(layerSymbols - decidedSymbols).^2);
        
        instantaneousSINR(nSlot+1, layerIdx) = 10*log10(signalPower / noisePower);
    end

    mesh(abs(estChGridLayers(:,:,1,1)));
    title('Channel Estimate');
    xlabel('OFDM Symbol');
    ylabel("Subcarrier");
    zlabel("Magnitude");

    constPlot.ChannelNames = "Layer "+(pdsch.NumLayers:-1:1);
    constPlot.ShowLegend = true;
    constPlot(fliplr(pdschEq));

    [dlschLLRs,rxSymbols] = nrPDSCHDecode(carrier,pdsch,pdschEq,noiseEst);

    csi = nrLayerDemap(csi);
    for cwIdx = 1:pdsch.NumCodewords
        Qm = length(dlschLLRs{cwIdx})/length(rxSymbols{cwIdx});
        csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);
        dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:);
    end

    decodeDLSCH.TransportBlockLength = trBlkSizes;
    [decbits,blkerr] = decodeDLSCH(dlschLLRs,pdsch.Modulation,pdsch.NumLayers, ...
        harqEntity.RedundancyVersion,harqEntity.HARQProcessID);

    if disableHARQ
        statusReport = "HARQ disabled - single transmission";
    else
        statusReport = updateAndAdvance(harqEntity,blkerr,trBlkSizes,pdschInfo.G);
    end

    disp("Slot "+(nSlot)+". "+statusReport);

    % ========== Track Performance Metrics ==========
    totalTransmissions = totalTransmissions + 1;

    for cwIdx = 1:pdsch.NumCodewords
        if disableHARQ || harqEntity.NewData(cwIdx)
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
            if disableHARQ || harqEntity.SequenceTimeout(cwIdx)
                blockErrors = blockErrors + 1;
            end
        end
    end

    perSlotBER(nSlot+1) = sum(blkerr) / pdsch.NumCodewords;

    disp("Slot "+(nSlot)+". "+statusReport);
end

% ========== Enhanced Results Display ==========
fprintf('\n========================================\n');
fprintf('     5G NR HARQ Simulation Results      \n');
fprintf('========================================\n');

fprintf('\n--- Configuration ---\n');
if disableHARQ
    fprintf('HARQ Type: Disabled\n');
    fprintf('RV Sequence: [0]\n');
    fprintf('Number of HARQ Processes: 1\n');
else
    fprintf('HARQ Type: Chase Combining\n');
    fprintf('RV Sequence: [%s]\n', num2str(rvSeq));
    fprintf('Number of HARQ Processes: %d\n', NHARQProcesses);
end
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

fprintf('\n--- Performance Metrics ---\n');
fprintf('Successful Blocks: %d / %d\n', successfulBlocks, totalInitialTransmissions);
fprintf('Failed Blocks (after max retx): %d\n', blockErrors);
fprintf('Block Error Rate (BLER): %.4f (%.2f%%)\n', ...
    blockErrors / totalInitialTransmissions, ...
    (blockErrors / totalInitialTransmissions) * 100);

fprintf('\n--- Throughput Analysis ---\n');
fprintf('Total Bits Transmitted: %d\n', totalTxBits);
fprintf('Total Bits Received (success): %d\n', totalRxBits);
fprintf('Throughput Efficiency: %.2f%%\n', (totalRxBits / totalTxBits) * 100);
fprintf('Effective Code Rate: %.4f\n', totalRxBits / (totalTransmissions * pdschInfo.G));

avgBitsPerSlot = totalRxBits / totalNoSlots;
slotDuration = carrier.SlotsPerSubframe / (carrier.SubcarrierSpacing / 15e3) * 1e-3;
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

coherenceBW_Hz = 1 / (5 * rmsDelaySpread);
fprintf('\nCoherence Bandwidth (50%% corr): %.2f kHz\n', coherenceBW_Hz / 1e3);

if isprop(channel, 'MaximumDopplerShift') && channel.MaximumDopplerShift > 0
    coherenceTime_s = 9 / (16 * pi * channel.MaximumDopplerShift);
    fprintf('Coherence Time (50%% corr): %.2f ms\n', coherenceTime_s * 1e3);
end

if logging == 1
    diary off;
end