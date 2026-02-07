% DL-SCH and PDSCH Transmit and Receive Processing Chain: 
% This example shows how to use 5G Toolbox™ features to model a 5G NR physical downlink shared channel (PDSCH) link, including all of the steps from transport block generation to bit decoding at the receiver end.

% OPTIONAL LOGGING VARIABLE: 
logging = 1;

% Test to Log the Command Window Output 
% logFile = sprintf('runlog_%s.txt', datestr(now,'yyyymmdd_HHMMSS'));
if logging == 1 
    logDir = 'tests/AustinS/logs';
    if ~exist(logDir,'dir'); mkdir(logDir); end
    logFile = fullfile(logDir, sprintf('runlog_%s.txt', datestr(now,'yyyymmdd_HHMMSS')));
    diary(logFile);
end

% SET TEST PARAMETERS HERE: 
% ------------------------------------------------------------------------------------------------
SNRdB = 8.0;
Modulation = "QPSK";                                                                               % must be 'QPSK', '16QAM', '64QAM', '256QAM', '1024QAM'
NHARQProcesses = 16;
rvSeq = [0 2 3 1];                                                                                  % originally set as [0 2 3 1]
% Specify the number of transmit and receive antennas 
nTxAnts = 8;
nRxAnts = 8;
% Specify Channel Model
DelayProfile = "TDL-C";

% Define a struct of parameters set in this simulation run 
% ------------------------------------------------------------------------------------------------
runParams.SNRdB = SNRdB;
runParams.Modulation = Modulation;
runParams.NHARQProcesses = NHARQProcesses;
runParams.rvSeq = rvSeq;
runParams.nTxAnts = nTxAnts;
runParams.nRxAnts = nRxAnts;
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
pdsch.NumLayers = 2;
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
fprintf('Git Commit: %s\n', strtrim(gitHash));
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

    % Estimate instantaneous SINR per layer from equalized symbols
    for layerIdx = 1:pdsch.NumLayers
        layerSymbols = pdschEq(:, layerIdx);
        refConstellation = getConstellationRefPoints(pdsch.Modulation);
        
        % Find nearest constellation point for each symbol
        [~, symbolDecisions] = min(abs(layerSymbols - refConstellation.'), [], 2);
        decidedSymbols = refConstellation(symbolDecisions);
        
        % Calculate signal and noise power
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
    % Constellation for the first layer has a higher SNR than that for the
    % last layer. Flip the layers so that the constellations do not mask
    % each other.
    constPlot(fliplr(pdschEq));

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

    disp("Slot "+(nSlot)+". "+statusReport);

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
end % for nSlot = 0:totalNoSlots

% % ========== Results Display ==========
% fprintf('\n========== Simulation Results ==========\n');
% fprintf('HARQ Configuration: Chase Combining (RV=[%s])\n', num2str(rvSeq));
% fprintf('SNR: %.2f dB\n', SNRdB);
% fprintf('Transport Blocks: %d\n', totalNoSlots);
% 
% % maxThroughput = sum(trBlkSizes);
% % totalNumRxBits = sum(decbits,2);
% 
% % fprintf('\n--- Overall Performance ---\n');
% % fprintf('Block Error Rate: %.4f (%.2f%%)\n', blkerr/totalNoSlots, ...
% %     blkerr/totalNoSlots*100);
% % fprintf('Throughput: %.2f%%\n', totalNumRxBits*100/maxThroughput);

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

% % ========== Performance Visualization ==========
% figure('Position', [100 100 1200 400]);
% 
% % Plot 1: Success rate over time
% subplot(1,3,1);
% plot(0:totalNoSlots-1, movmean(perSlotSuccess, 5), 'LineWidth', 2);
% grid on;
% xlabel('Slot Number');
% ylabel('Success Rate (5-slot moving avg)');
% title('Decoding Success Rate vs. Time');
% ylim([0 1.1]);
% 
% % Plot 2: Cumulative throughput
% subplot(1,3,2);
% cumulativeBits = cumsum(perSlotSuccess .* trBlkSizes(1));
% plot(0:totalNoSlots-1, cumulativeBits / 1e6, 'LineWidth', 2);
% grid on;
% xlabel('Slot Number');
% ylabel('Throughput (Mbits)');
% title('Cumulative Throughput');
% 
% % Plot 3: HARQ retransmission distribution
% subplot(1,3,3);
% bar([totalInitialTransmissions, totalRetransmissions]);
% set(gca, 'XTickLabel', {'Initial TX', 'Retransmissions'});
% ylabel('Count');
% title('Transmission Distribution');
% grid on;


% 
% % ========== Channel Performance Visualization ==========
% figure('Position', [100 100 1400 800]);
% 
% % Plot 1: Channel gain over time
% subplot(2,3,1);
% plot(0:totalNoSlots-1, avgChannelGain, 'LineWidth', 1.5);
% grid on;
% xlabel('Slot Number');
% ylabel('Channel Gain (dB)');
% title('Average Channel Gain vs. Time');
% 
% % Plot 2: SINR per layer
% subplot(2,3,2);
% hold on;
% for layerIdx = 1:pdsch.NumLayers
%     plot(0:totalNoSlots-1, instantaneousSINR(:, layerIdx), ...
%         'LineWidth', 1.5, 'DisplayName', sprintf('Layer %d', layerIdx));
% end
% hold off;
% grid on;
% xlabel('Slot Number');
% ylabel('SINR (dB)');
% title('Instantaneous SINR per Layer');
% legend('show');
% 
% % Plot 3: Condition number
% subplot(2,3,3);
% plot(0:totalNoSlots-1, conditionNumber, 'LineWidth', 1.5);
% grid on;
% xlabel('Slot Number');
% ylabel('Condition Number');
% title('MIMO Channel Condition Number');
% 
% % Plot 4: Power delay profile
% subplot(2,3,4);
% stem(chInfo.PathDelays * 1e9, chInfo.AveragePathGains, 'LineWidth', 2, 'MarkerSize', 8);
% grid on;
% xlabel('Delay (ns)');
% ylabel('Average Path Gain (dB)');
% title(sprintf('%s Power Delay Profile', channel.DelayProfile));
% 
% % Plot 5: SINR histogram
% subplot(2,3,5);
% histogram(instantaneousSINR(:), 20, 'FaceAlpha', 0.7);
% grid on;
% xlabel('SINR (dB)');
% ylabel('Frequency');
% title('SINR Distribution (All Layers)');
% 
% % Plot 6: Channel gain vs Block errors
% subplot(2,3,6);
% scatter(avgChannelGain, perSlotBER, 50, 'filled');
% grid on;
% xlabel('Channel Gain (dB)');
% ylabel('Block Error Rate');
% title('Channel Gain vs. Block Error Rate');



% Final line for Command Window tracking
if logging == 1
    diary off;
end