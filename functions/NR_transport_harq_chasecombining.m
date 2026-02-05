% 5G NR Link-Level Simulation with Chase Combining
% Using MATLAB's nrDLSCHDecoder Native Chase Combining Support
%
% KEY: Just use RV=[0 0 0 0] and let the decoder handle everything!
% The decoder automatically does soft combining when it sees repeated RVs.
% No manual LLR buffer management needed.

% Step 1: 
noTransportBlocks = 100;
SNRdB = 0.5; % SNR in dB

% Step 2: 
rng(42);

% Step 3: 
% DL-SCH parameters
codeRate = 490/1024;
NHARQProcesses = 16;

% ============= KEY DIFFERENCE =============
% Chase Combining: Use RV=0 for all retransmissions
% The decoder will automatically perform soft combining
rvSeq = [0 0 0 0];  % Chase Combining

% For comparison, Incremental Redundancy would use:
% rvSeq = [0 2 3 1];  % Incremental Redundancy
% ==========================================

% Step 4: 
% Create DL-SCH encoder and decoder objects
encodeDLSCH = nrDLSCH;
encodeDLSCH.MultipleHARQProcesses = true;
encodeDLSCH.TargetCodeRate = codeRate;

decodeDLSCH = nrDLSCHDecoder;
decodeDLSCH.MultipleHARQProcesses = true;
decodeDLSCH.TargetCodeRate = codeRate;
decodeDLSCH.LDPCDecodingAlgorithm = "Normalized min-sum";
decodeDLSCH.MaximumLDPCIterationCount = 6;

% Step 5: Carrier configuration
SCS = 15;
NRB = 52;

carrier = nrCarrierConfig;
carrier.NSizeGrid = NRB;
carrier.SubcarrierSpacing = SCS;
carrier.CyclicPrefix = "Normal";

% Step 6: PDSCH configuration
modulation = "16QAM";

pdsch = nrPDSCHConfig;
pdsch.Modulation = modulation;
pdsch.PRBSet = 0:NRB-1;
pdsch.NumLayers = 1;

% Step 7: HARQ Management
harqEntity = HARQEntity(0:NHARQProcesses-1,rvSeq,pdsch.NumCodewords);

% Step 7b: Statistics tracking for Chase Combining analysis
ccStats = struct();
ccStats.transmissionsPerTB = [];
ccStats.firstTxSuccess = 0;
ccStats.successAfterRetx = 0;
ccStats.totalRetransmissions = 0;

% Step 8: BER Simulation
noiseVar = 1./(10.^(SNRdB/10));
numBlkErr = 0;
numRxBits = [];
txedTrBlkSizes = [];

fprintf('========== Starting Simulation ==========\n');
fprintf('Configuration:\n');
fprintf('  HARQ Type: Chase Combining (RV sequence: [%s])\n', num2str(rvSeq));
fprintf('  SNR: %.2f dB\n', SNRdB);
fprintf('  Transport Blocks: %d\n', noTransportBlocks);
fprintf('  HARQ Processes: %d\n', NHARQProcesses);
fprintf('=========================================\n\n');

% Track current transmission count per process
currentTxCount = zeros(NHARQProcesses, 1);

for nTrBlk = 1:noTransportBlocks
    carrier.NSlot = carrier.NSlot+1;

    % Calculate transport block size
    [~,pdschInfo] = nrPDSCHIndices(carrier,pdsch);
    Xoh_PDSCH = 0;
    trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet),...
                       pdschInfo.NREPerPRB,codeRate,Xoh_PDSCH);
 
    % HARQ Processing
    for cwIdx = 1:pdsch.NumCodewords
        if harqEntity.NewData(cwIdx)
            % New transport block
            trBlk = randi([0 1],trBlkSizes(cwIdx),1);
            setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);

            if harqEntity.SequenceTimeout(cwIdx)
                % Previous TB failed after all retransmissions
                resetSoftBuffer(decodeDLSCH,cwIdx-1,harqEntity.HARQProcessID);
                
                % Record that this TB needed max transmissions
                ccStats.transmissionsPerTB = [ccStats.transmissionsPerTB length(rvSeq)];
                currentTxCount(harqEntity.HARQProcessID+1) = 0;
            else
                % Reset counter for new TB
                currentTxCount(harqEntity.HARQProcessID+1) = 0;
            end
        end
    end
    
    % Increment transmission counter for current process
    currentProcessID = harqEntity.HARQProcessID + 1;
    currentTxCount(currentProcessID) = currentTxCount(currentProcessID) + 1;
    isRetransmission = currentTxCount(currentProcessID) > 1;

    % DL-SCH Encoding
    % The encoder generates the same coded bits for each RV=0 transmission
    codedTrBlock = encodeDLSCH(pdsch.Modulation,pdsch.NumLayers,pdschInfo.G, ...
        harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
    
    % PDSCH Encoding
    modOut = nrPDSCH(carrier,pdsch,codedTrBlock);
    
    % AWGN Channel
    rxSig = awgn(modOut,SNRdB);    

    % PDSCH Demodulation
    rxLLR = nrPDSCHDecode(carrier,pdsch,rxSig,noiseVar);

    % ============= DL-SCH DECODING =============
    % The decoder automatically performs Chase Combining:
    % - On first transmission (RV=0): Stores LLRs in soft buffer
    % - On retransmissions (RV=0): Adds new LLRs to soft buffer
    % - This is all handled internally by the decoder!
    decodeDLSCH.TransportBlockLength = trBlkSizes;
    [decbits,blkerr] = decodeDLSCH(rxLLR,pdsch.Modulation,pdsch.NumLayers, ...
        harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
    % ===========================================
    
    % Track statistics
    if ~blkerr
        % Successful decode
        numTx = currentTxCount(currentProcessID);
        ccStats.transmissionsPerTB = [ccStats.transmissionsPerTB numTx];
        
        if numTx == 1
            ccStats.firstTxSuccess = ccStats.firstTxSuccess + 1;
        else
            ccStats.successAfterRetx = ccStats.successAfterRetx + 1;
            ccStats.totalRetransmissions = ccStats.totalRetransmissions + (numTx - 1);
        end
        
        % Reset counter on success
        currentTxCount(currentProcessID) = 0;
    else
        % Failed decode
        if isRetransmission
            ccStats.totalRetransmissions = ccStats.totalRetransmissions + 1;
        end
    end
    
    % RESULTS
    if(any(trBlkSizes ~= 0))
        numRxBits = [numRxBits trBlkSizes.*(1-blkerr)];
        txedTrBlkSizes = [txedTrBlkSizes trBlkSizes];
    end
    
    if blkerr   
        numBlkErr = numBlkErr + 1;
    end

    % HARQ Update
    statusReport = updateAndAdvance(harqEntity,blkerr,trBlkSizes,pdschInfo.G);    

    % Display information
    txNum = currentTxCount(currentProcessID);
    if txNum > 1
        theoreticalGain = 10*log10(txNum);
        ccInfo = sprintf(" [Chase: Tx#%d, ~%.1f dB gain]", txNum, theoreticalGain);
    else
        ccInfo = " [Chase: First Tx]";
    end
    
    % Show periodic updates
    if mod(nTrBlk, 10) == 0 || nTrBlk <= 5
        disp("Slot "+(nTrBlk)+". "+statusReport+ccInfo);
    end
   
end

% ========== Results Display ==========
fprintf('\n========== Simulation Results ==========\n');
fprintf('HARQ Configuration: Chase Combining (RV=[%s])\n', num2str(rvSeq));
fprintf('SNR: %.2f dB\n', SNRdB);
fprintf('Transport Blocks: %d\n', noTransportBlocks);

maxThroughput = sum(txedTrBlkSizes);
totalNumRxBits = sum(numRxBits,2);

fprintf('\n--- Overall Performance ---\n');
fprintf('Block Error Rate: %.4f (%.2f%%)\n', numBlkErr/noTransportBlocks, ...
    numBlkErr/noTransportBlocks*100);
fprintf('Throughput: %.2f%%\n', totalNumRxBits*100/maxThroughput);

fprintf('\n--- Chase Combining Statistics ---\n');
fprintf('Success on First Transmission: %d (%.1f%%)\n', ...
    ccStats.firstTxSuccess, ccStats.firstTxSuccess/noTransportBlocks*100);
fprintf('Success After Retransmission(s): %d (%.1f%%)\n', ...
    ccStats.successAfterRetx, ccStats.successAfterRetx/noTransportBlocks*100);
fprintf('Total Failed (after max retx): %d (%.1f%%)\n', ...
    numBlkErr, numBlkErr/noTransportBlocks*100);
fprintf('Total Retransmissions: %d\n', ccStats.totalRetransmissions);

if ccStats.successAfterRetx > 0
    fprintf('Average Retransmissions (when needed): %.2f\n', ...
        ccStats.totalRetransmissions/ccStats.successAfterRetx);
end

% Transmission distribution
if ~isempty(ccStats.transmissionsPerTB)
    fprintf('\n--- Transmission Distribution ---\n');
    fprintf('Average Transmissions per Successful TB: %.2f\n', ...
        mean(ccStats.transmissionsPerTB));
    
    for tx = 1:max(ccStats.transmissionsPerTB)
        count = sum(ccStats.transmissionsPerTB == tx);
        if count > 0
            pct = count/length(ccStats.transmissionsPerTB)*100;
            gainDB = 10*log10(tx);
            fprintf('  %d transmission(s): %d TBs (%.1f%%) [~%.1f dB effective gain]\n', ...
                tx, count, pct, gainDB);
        end
    end
end

% Theory vs Practice
fprintf('\n--- Theoretical Analysis ---\n');
fprintf('Chase Combining provides ~3 dB gain per retransmission:\n');
fprintf('  Base SNR: %.2f dB\n', SNRdB);
for n = 1:4
    effSNR = SNRdB + 10*log10(n);
    fprintf('  After %d transmissions: Effective SNR ≈ %.2f dB (gain: %.1f dB)\n', ...
        n, effSNR, 10*log10(n));
end

fprintf('\n========================================\n');

fprintf('\nNOTE: The nrDLSCHDecoder automatically performs Chase Combining\n');
fprintf('when the same RV is used for retransmissions. The soft buffer\n');
fprintf('management and LLR combining are handled internally.\n');