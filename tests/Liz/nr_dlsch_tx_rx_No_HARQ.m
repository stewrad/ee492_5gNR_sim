% 5G NR PDSCH Link Simulation - NO HARQ
% Outputs only SNR and BLER for charting

% --------------------------------------------------------------------------------
% USER SETTINGS
% --------------------------------------------------------------------------------
SNRdBList = 0:2:12;                  % SNR sweep for chart
Modulation = "256QAM";               % 'QPSK', '16QAM', '64QAM', '256QAM', '1024QAM'
nTxAnts = 8;
nRxAnts = 8;
NumLayers = 2;                       % keep <= min(nTxAnts,nRxAnts)
DelayProfile = "TDL-C";
totalNoSlots = 50;                  % increase for smoother BLER curve
perfectEstimation = false;
rng("default");

% --------------------------------------------------------------------------------
% CARRIER CONFIGURATION
% --------------------------------------------------------------------------------
carrier = nrCarrierConfig;

pdsch = nrPDSCHConfig;
pdsch.Modulation = Modulation;
pdsch.NumLayers = NumLayers;
pdsch.PRBSet = 0:carrier.NSizeGrid-1;

pdsch.DMRS.DMRSAdditionalPosition = 1;
pdsch.DMRS.DMRSConfigurationType = 1;
pdsch.DMRS.DMRSLength = 2;

if pdsch.NumLayers > min(nTxAnts,nRxAnts)
    error("The number of layers (" + string(pdsch.NumLayers) + ...
          ") must be <= min(nTxAnts,nRxAnts) (" + string(min(nTxAnts,nRxAnts)) + ").");
end

% For NumLayers <= 4, there is 1 codeword
if pdsch.NumCodewords ~= 1
    error("This simplified no-HARQ script assumes one codeword. Reduce NumLayers if needed.");
end

codeRate = 490/1024;
rv = 0;   % no HARQ, single transmission only

% --------------------------------------------------------------------------------
% CHANNEL CONFIGURATION
% --------------------------------------------------------------------------------
channel = nrTDLChannel;
channel.DelayProfile = DelayProfile;
channel.NumTransmitAntennas = nTxAnts;
channel.NumReceiveAntennas = nRxAnts;

ofdmInfo = nrOFDMInfo(carrier);
channel.SampleRate = ofdmInfo.SampleRate;
channel.ChannelResponseOutput = 'ofdm-response';

% Initial precoder estimate
estChannelGrid = getInitialChannelEstimate(channel,carrier);
newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChannelGrid);

% --------------------------------------------------------------------------------
% DL-SCH ENCODER / DECODER (NO HARQ)
% --------------------------------------------------------------------------------
encodeDLSCH = nrDLSCH;
encodeDLSCH.MultipleHARQProcesses = false;
encodeDLSCH.TargetCodeRate = codeRate;

decodeDLSCH = nrDLSCHDecoder;
decodeDLSCH.MultipleHARQProcesses = false;
decodeDLSCH.TargetCodeRate = codeRate;
decodeDLSCH.LDPCDecodingAlgorithm = "Normalized min-sum";
decodeDLSCH.MaximumLDPCIterationCount = 6;

% --------------------------------------------------------------------------------
% RESULTS STORAGE
% --------------------------------------------------------------------------------
results = zeros(numel(SNRdBList),2);   % col 1 = SNR, col 2 = BLER

% --------------------------------------------------------------------------------
% SNR SWEEP
% --------------------------------------------------------------------------------
for snrIdx = 1:numel(SNRdBList)

    SNRdB = SNRdBList(snrIdx);

    failedBlocks = 0;
    totalBlocks = 0;

    reset(channel);

    for nSlot = 0:totalNoSlots-1
        carrier.NSlot = nSlot;

        [pdschIndices,pdschInfo] = nrPDSCHIndices(carrier,pdsch);

        Xoh_PDSCH = 0;
        trBlkSize = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet), ...
                          pdschInfo.NREPerPRB,codeRate,Xoh_PDSCH);

        % Fresh block every slot, no HARQ memory
        trBlk = randi([0 1],trBlkSize,1);
        setTransportBlock(encodeDLSCH,trBlk,0);

        % Reset decoder soft buffer every slot for true no-HARQ behavior
        resetSoftBuffer(decodeDLSCH,0);

        codedTrBlock = encodeDLSCH(pdsch.Modulation,pdsch.NumLayers,pdschInfo.G,rv);

        pdschSymbols = nrPDSCH(carrier,pdsch,codedTrBlock);

        precodingWeights = newPrecodingWeight;
        pdschSymbolsPrecoded = pdschSymbols * precodingWeights;

        dmrsSymbols = nrPDSCHDMRS(carrier,pdsch);
        dmrsIndices = nrPDSCHDMRSIndices(carrier,pdsch);

        pdschGrid = nrResourceGrid(carrier,nTxAnts);

        [~,pdschAntIndices] = nrExtractResources(pdschIndices,pdschGrid);
        pdschGrid(pdschAntIndices) = pdschSymbolsPrecoded;

        for p = 1:size(dmrsSymbols,2)
            [~,dmrsAntIndices] = nrExtractResources(dmrsIndices(:,p),pdschGrid);
            pdschGrid(dmrsAntIndices) = pdschGrid(dmrsAntIndices) + ...
                                        dmrsSymbols(:,p) * precodingWeights(p,:);
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
            [offset,~] = nrTimingEstimate(carrier,rxWaveform,dmrsIndices,dmrsSymbols);
        end

        offset = max(0,offset);
        rxWaveform = rxWaveform(1+offset:end,:);

        rxGrid = nrOFDMDemodulate(carrier,rxWaveform);

        if perfectEstimation
            estChGridAnts = ofdmChannelResponse;
            noiseEst = nVar;
            newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChGridAnts);
            estChGridLayers = precodeChannelEstimate(estChGridAnts,precodingWeights.');
        else
            [estChGridLayers,noiseEst] = nrChannelEstimate( ...
                carrier,rxGrid,dmrsIndices,dmrsSymbols, ...
                'CDMLengths',pdsch.DMRS.CDMLengths);

            estChGridAnts = precodeChannelEstimate(estChGridLayers,conj(precodingWeights));
            newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChGridAnts);
        end

        [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxGrid,estChGridLayers);
        [pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);

        [dlschLLRs,rxSymbols] = nrPDSCHDecode(carrier,pdsch,pdschEq,noiseEst);

        csi = nrLayerDemap(csi);
        Qm = length(dlschLLRs{1}) / length(rxSymbols{1});
        csi{1} = repmat(csi{1}.',Qm,1);
        dlschLLRs{1} = dlschLLRs{1} .* csi{1}(:);

        decodeDLSCH.TransportBlockLength = trBlkSize;
        [~,blkerr] = decodeDLSCH(dlschLLRs,pdsch.Modulation,pdsch.NumLayers,rv);

        failedBlocks = failedBlocks + blkerr;
        totalBlocks = totalBlocks + 1;
    end

    BLER = failedBlocks / totalBlocks;

    results(snrIdx,1) = SNRdB;
    results(snrIdx,2) = BLER;

    fprintf('SNR = %2.1f dB, BLER = %.6f\n', SNRdB, BLER);
end

% --------------------------------------------------------------------------------
% FINAL OUTPUT TABLE
% --------------------------------------------------------------------------------
resultsTable = array2table(results, 'VariableNames', {'SNRdB','BLER'});
disp(resultsTable);