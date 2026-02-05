% Step 1: 
% Specify the number of transport blocks and the SNR 
noTransportBlocks = 100;
SNRdB = 1.5; % SNR in dB

% Step 2: 
% Reset random number generator for reproducibility
rng(42);

% Step 3: 
% DL-SCH parameters
% Specify the code rate, the number of HARQ processes, and the redundancy values (RVs) sequence. This sequence controls the redundancy version retransmissions in case of error. 
codeRate = 490/1024;
NHARQProcesses = 16; % Number of parallel HARQ processes to use
rvSeq = [0 2 3 1];

% Step 4: 
% Create the DL-SCH encoder and decoder objects. To use multiple processes, set the MultipleHARQProcesses property to true for both objects. To enable retransmissions for multiple HARQ processes, the encoder buffers the input bits. The decoder needs a similar mechanism to enable soft combining of retransmissions for each HARQ process. 
% Create DL-SCH encoder object
encodeDLSCH = nrDLSCH;
encodeDLSCH.MultipleHARQProcesses = true;
encodeDLSCH.TargetCodeRate = codeRate;
% Create DL-SCH decoder object
decodeDLSCH = nrDLSCHDecoder;
decodeDLSCH.MultipleHARQProcesses = true;
decodeDLSCH.TargetCodeRate = codeRate;
decodeDLSCH.LDPCDecodingAlgorithm = "Normalized min-sum";
decodeDLSCH.MaximumLDPCIterationCount = 6;
% % The DL-SCH encoder and decoder objects can model up to 16 HARQ processes. The encoder and decoder objects use the HARQprocessID property of the HARQ entity object to identify the active HARQ process when performing any of these operations.
% Setting new transport block to transmit
% Encoding data
% Resetting soft buffers
% Decoding data

% Step 5: 
% Specify the carrier and PDSCH parameters. These parameters are used for PDSCH encoding and decoding and for calculating the transport block size.
% Numerology
SCS = 15;                         % SCS: 15, 30, 60, 120 or 240 (kHz)
NRB = 52;                         % BW in number of RBs (52 RBs at 15 kHz SCS for 10 MHz BW)
% Create a carrier object, specifying the subcarrier spacing (SCS) and the bandwidth (BW).
carrier = nrCarrierConfig;
carrier.NSizeGrid = NRB;
carrier.SubcarrierSpacing = SCS;
carrier.CyclicPrefix = "Normal";  % "Normal" or "Extended"

% Step 6: 
% Create a  PDSCH configuration object. The PDSCH parameters determine the available bit capacity and the transport block size.
modulation = "16QAM";             % Modulation scheme

pdsch = nrPDSCHConfig;
pdsch.Modulation = modulation;
pdsch.PRBSet = 0:NRB-1;           % Assume full band allocation
pdsch.NumLayers = 1;              % Assume only one layer and one codeword


% Step 7: 
% HARQ Management ------------------
% Create a HARQ entity object to manage the HARQ processes. For each HARQ processes, the object stores these elements:
% HARQ ID number.
% RV.
% Transmission number, which indicates how many times a certain transport block has been transmitted.
% Flag to indicate whether new data is required. New data is required when a transport block is received successfully or if a sequence timeout has occurred (all RV transmissions have failed).
% Flag to indicate whether a sequence timeout has occurred (all RV transmissions have failed).
harqEntity = HARQEntity(0:NHARQProcesses-1,rvSeq,pdsch.NumCodewords);
% The HARQ entity is used to manage the buffers in the DL-SCH encoder and decoder.

% Step 8: 
% BER Simulation -------------------
% Loop over a number of transport blocks. For each transport block:
% Calculate the transport block size in number of bits.
% Generate new data block or reset buffers in the decoder.
% Apply DL-SCH encoding.
% Modulate bits to symbols.
% Apply AWGN.
% Demodulate soft bits (symbols to soft bits).
% Decode the DL-SCH.
% Update the HARQ processes.
% Initialize loop variables
noiseVar = 1./(10.^(SNRdB/10)); % Noise variance
numBlkErr = 0;                  % Number of block errors
numRxBits = [];                 % Number of successfully received bits per transmission
txedTrBlkSizes = [];            % Number of transmitted info bits per transmission
for nTrBlk = 1:noTransportBlocks
    % A transport block or transmission time interval (TTI) corresponds to
    % one slot
    carrier.NSlot = carrier.NSlot+1;

    % Calculate the transport block size. -------------------------------
    % Generate PDSCH indices info, which is used to calculate the transport
    % block size
    [~,pdschInfo] = nrPDSCHIndices(carrier,pdsch);
    % Calculate transport block sizes
    Xoh_PDSCH = 0;
    trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet),pdschInfo.NREPerPRB,codeRate,Xoh_PDSCH);
 
    % HARQ Processing (Buffer Management)
    % This section explains the buffer management in the encoder and decoder.
    % DL-SCH encoder buffers: Generate a new transport block if new data is required for the active HARQ process. Store the transport block in the corresponding buffer. If no new data is required, the buffered bits in the DL-SCH encoder are used for retransmission.
    % DL-SCH decoder buffers: The soft buffers in the receiver store previously received versions of the same transport block. These buffers are cleared automatically upon successful reception (no CRC error). However, if the RV sequence ends without successful decoding, the buffers must be flushed manually by calling the resetSoftBuffer object function.
    % Get new transport blocks and flush decoder soft buffer, as required
    for cwIdx = 1:pdsch.NumCodewords
        if harqEntity.NewData(cwIdx)
            % Create and store a new transport block for transmission
            trBlk = randi([0 1],trBlkSizes(cwIdx),1);
            setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);

            % If the previous RV sequence ends without successful decoding,
            % flush the soft buffer explicitly
            if harqEntity.SequenceTimeout(cwIdx)
                resetSoftBuffer(decodeDLSCH,cwIdx-1,harqEntity.HARQProcessID);
            end
        end
    end

    % DL-SCH Encoding: Encode DL-SCH transport blocks
    codedTrBlock = encodeDLSCH(pdsch.Modulation,pdsch.NumLayers,pdschInfo.G, ...
        harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
    % PDSCH Encoding: Generate PDSCH symbols
    modOut = nrPDSCH(carrier,pdsch,codedTrBlock);
    
    % AWGN Channel
    rxSig = awgn(modOut,SNRdB);    

   % PDSCH Demodulation: 
   % Soft demodulate the received symbols
    rxLLR = nrPDSCHDecode(carrier,pdsch,rxSig,noiseVar);

    % DL-SCH Decoding: Apply DL-SCH decoding
    decodeDLSCH.TransportBlockLength = trBlkSizes;
    [decbits,blkerr] = decodeDLSCH(rxLLR,pdsch.Modulation,pdsch.NumLayers, ...
        harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
    
    % RESULTS ----
    % Store values to calculate throughput (only for active transport blocks)
    if(any(trBlkSizes ~= 0))
        numRxBits = [numRxBits trBlkSizes.*(1-blkerr)];
        txedTrBlkSizes = [txedTrBlkSizes trBlkSizes];
    end
    
    if blkerr   
        numBlkErr = numBlkErr + 1;
    end

    % HARQ Process Update 
    % Update the current HARQ process with the CRC error, and then advance to the next process. This step updates the information related to the active HARQ process in the HARQ entity.
    statusReport = updateAndAdvance(harqEntity,blkerr,trBlkSizes,pdschInfo.G);    

    % Display information about current decoding attempt
    disp("Slot "+(nTrBlk)+". "+statusReport);
   
end % for nTrBlk = 1:noTransportBlocks


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