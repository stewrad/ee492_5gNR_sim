SNRvec = 0:2:12;
finalBLER_vs_SNR = zeros(size(SNRvec));

for iSNR = 1:length(SNRvec)

    SNRdB = SNRvec(iSNR);

    % reset per-SNR counters
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

    instantaneousSINR = zeros(totalNoSlots, pdsch.NumLayers);
    avgChannelGain = zeros(totalNoSlots, 1);
    conditionNumber = zeros(totalNoSlots, 1);

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

    reset(decodeDLSCH);

    harqEntity = HARQEntity(0:NHARQProcesses-1,rvSeq,pdsch.NumCodewords);

    for nSlot = 0:totalNoSlots-1
        carrier.NSlot = nSlot;

        [pdschIndices,pdschInfo] = nrPDSCHIndices(carrier,pdsch);

        Xoh_PDSCH = 0;
        trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet),pdschInfo.NREPerPRB,codeRate,Xoh_PDSCH);

        currentIsNewData = false(1,pdsch.NumCodewords);
        currentTimedOut  = false(1,pdsch.NumCodewords);

        for cwIdx = 1:pdsch.NumCodewords
            currentIsNewData(cwIdx) = harqEntity.NewData(cwIdx);
            currentTimedOut(cwIdx)  = harqEntity.SequenceTimeout(cwIdx);

            if currentIsNewData(cwIdx)
                trBlk = randi([0 1],trBlkSizes(cwIdx),1);
                setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);
                totalInitialTransmissions = totalInitialTransmissions + 1;

                if currentTimedOut(cwIdx)
                    resetSoftBuffer(decodeDLSCH,cwIdx-1,harqEntity.HARQProcessID);
                end
            else
                totalRetransmissions = totalRetransmissions + 1;
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
        [decbits,blkerr] = decodeDLSCH(dlschLLRs,pdsch.Modulation,pdsch.NumLayers,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);

        totalTransmissions = totalTransmissions + pdsch.NumCodewords;
        perfState.attemptsTotal = perfState.attemptsTotal + pdsch.NumCodewords;
        perfState.attemptsFailed = perfState.attemptsFailed + sum(blkerr);

        for cwIdx = 1:pdsch.NumCodewords
            if currentIsNewData(cwIdx) && ~blkerr(cwIdx)
                successfulBlocks = successfulBlocks + 1;
            end
            if currentTimedOut(cwIdx) && blkerr(cwIdx)
                blockErrors = blockErrors + 1;
            end
        end

        harqEntity.updateAndAdvance(blkerr,trBlkSizes,pdschInfo.G);
    end

    finalBLER_vs_SNR(iSNR) = blockErrors / max(1,totalInitialTransmissions);
    fprintf('SNR = %2d dB, Final BLER = %.4f\n', SNRdB, finalBLER_vs_SNR(iSNR));

end

figure;
plot(SNRvec, finalBLER_vs_SNR*100, '-o', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('BLER (%)');
title('BLER vs SNR');