%% matlabtrial5G_HARQ_selfcontained.m
% Self-contained 5G NR DL-SCH + PDSCH link simulation WITH HARQ
%
% Implements:
% - 16 parallel HARQ processes (IDs 0..15)
% - RV sequence [0 2 3 1]
% - New data on ACK or RV sequence timeout
% - Soft combining across retransmissions via nrDLSCHDecoder buffers
% - Soft buffer flush on timeout  
%
% Includes:
% - 2-layer 16QAM PDSCH over 8x8 TDL-C channel
% - DM-RS-based timing and channel estimation (practical), MMSE equalization
% - SVD-based precoding updated slot-to-slot
%
% Requirements:
% - 5G Toolbox
% - Communications Toolbox (for comm.ConstellationDiagram)

clear; clc;

%% Simulation parameters
SNRdB = 10;                % SNR in dB (per RE per RX antenna)
totalNoSlots = 20;         % Number of slots to simulate
perfectEstimation = false; % Perfect sync + channel estimation flag
rng("default");            % Repeatability

%% Carrier configuration
carrier = nrCarrierConfig; % default numerology and grid
disp(carrier);

%% PDSCH and DM-RS configuration
pdsch = nrPDSCHConfig;
pdsch.Modulation = "16QAM";
pdsch.NumLayers  = 2;
pdsch.PRBSet     = 0:carrier.NSizeGrid-1; % full-band allocation

pdsch.DMRS.DMRSAdditionalPosition  = 1;
pdsch.DMRS.DMRSConfigurationType   = 1;
pdsch.DMRS.DMRSLength              = 2;
disp(pdsch.DMRS);

%% DL-SCH configuration
NHARQProcesses = 16;   % up to 16 supported by nrDLSCH/nrDLSCHDecoder
rvSeq = [0 2 3 1];

if pdsch.NumCodewords == 1
    codeRate = 490/1024;
else
    codeRate = [490 490]./1024;
end

% DL-SCH encoder
encodeDLSCH = nrDLSCH;
encodeDLSCH.MultipleHARQProcesses = true;
encodeDLSCH.TargetCodeRate = codeRate;

% DL-SCH decoder
decodeDLSCH = nrDLSCHDecoder;
decodeDLSCH.MultipleHARQProcesses = true;
decodeDLSCH.TargetCodeRate = codeRate;
decodeDLSCH.LDPCDecodingAlgorithm = "Normalized min-sum";
decodeDLSCH.MaximumLDPCIterationCount = 6;

%% HARQ manager (self-contained replacement for HARQEntity)
% Note: MATLAB arrays are 1-based. We store HARQ state in rows 1..NHARQProcesses,
% while the Toolbox HARQ process IDs are 0..NHARQProcesses-1.
harq = initHARQManager(NHARQProcesses, rvSeq, pdsch.NumCodewords);

%% Channel configuration
nTxAnts = 8;
nRxAnts = 8;

if pdsch.NumLayers > min(nTxAnts,nRxAnts)
    error("The number of layers must be <= min(nTxAnts,nRxAnts).");
end

channel = nrTDLChannel;
channel.DelayProfile = "TDL-C";
channel.NumTransmitAntennas = nTxAnts;
channel.NumReceiveAntennas  = nRxAnts;

ofdmInfo = nrOFDMInfo(carrier);
channel.SampleRate = ofdmInfo.SampleRate;

% Needed to support perfect channel response output and timing offset output
channel.ChannelResponseOutput = "ofdm-response";

%% Constellation diagram
constPlot = comm.ConstellationDiagram;
constPlot.ReferenceConstellation = getConstellationRefPoints(pdsch.Modulation);
constPlot.EnableMeasurements = true;

%% Initial timing offset and initial precoder
offset = 0;

% Initial (perfect) channel estimate for first-slot precoder
estChGrid0 = getInitialChannelEstimate(channel, carrier);
newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet, pdsch.NumLayers, estChGrid0);

%% Metrics
slotAnyBlkErr = false(totalNoSlots,1);

%% Slot loop
for nSlot = 0:totalNoSlots-1
    carrier.NSlot = nSlot;

    % ---- Select HARQ process (round-robin) ----
    pid0 = mod(nSlot, NHARQProcesses); % 0-based HARQ process ID for toolbox calls
    pidx = pid0 + 1;                  % 1-based index for harq state arrays

    % ---- PDSCH indices and transport block size ----
    [pdschIndices,pdschInfo] = nrPDSCHIndices(carrier, pdsch);

    Xoh_PDSCH = 0;
    trBlkSizes = nrTBS(pdsch.Modulation, pdsch.NumLayers, numel(pdsch.PRBSet), ...
        pdschInfo.NREPerPRB, codeRate, Xoh_PDSCH);

    % ---- HARQ buffer management (per codeword) ----
    for cw = 1:pdsch.NumCodewords
        if harq.NewData(pidx,cw)
            % Create/store new TB in encoder buffer for this (pid,cw)
            trBlk = randi([0 1], trBlkSizes(cw), 1);
            setTransportBlock(encodeDLSCH, trBlk, cw-1, pid0);

            % If previous RV sequence timed out, flush decoder soft buffer
            if harq.SequenceTimeout(pidx,cw)
                resetSoftBuffer(decodeDLSCH, cw-1, pid0);
                harq.SequenceTimeout(pidx,cw) = false;
            end
        end
    end

    % Current redundancy versions (vector length = NumCodewords)
    rv = harq.RedundancyVersion(pidx,:);

    % ---- DL-SCH encode (LDPC + rate match) ----
    codedTrBlock = encodeDLSCH(pdsch.Modulation, pdsch.NumLayers, pdschInfo.G, rv, pid0);

    % ---- PDSCH modulation ----
    pdschSymbols = nrPDSCH(carrier, pdsch, codedTrBlock);

    % ---- Precoding ----
    precodingWeights = newPrecodingWeight;          % NLayers x nTxAnts
    pdschSymbolsPrecoded = pdschSymbols * precodingWeights;

    % ---- DM-RS generation ----
    dmrsSymbols = nrPDSCHDMRS(carrier, pdsch);
    dmrsIndices = nrPDSCHDMRSIndices(carrier, pdsch);

    % ---- Map to resource grid (antenna-domain) ----
    pdschGrid = nrResourceGrid(carrier, nTxAnts);

    % Convert layer indices -> antenna indices since symbols already precoded
    [~,pdschAntIndices] = nrExtractResources(pdschIndices, pdschGrid);
    pdschGrid(pdschAntIndices) = pdschSymbolsPrecoded;

    % DM-RS precoding + mapping
    for p = 1:size(dmrsSymbols,2)
        [~,dmrsAntIndices] = nrExtractResources(dmrsIndices(:,p), pdschGrid);
        pdschGrid(dmrsAntIndices) = pdschGrid(dmrsAntIndices) + dmrsSymbols(:,p) * precodingWeights(p,:);
    end

    % ---- OFDM modulation ----
    [txWaveform,waveformInfo] = nrOFDMModulate(carrier, pdschGrid);

    % ---- Channel padding for delay ----
    chInfo = info(channel);
    maxChDelay = chInfo.MaximumChannelDelay;
    txWaveform = [txWaveform; zeros(maxChDelay, size(txWaveform,2))];

    % ---- Propagation channel + noise ----
    [rxWaveform, ofdmChannelResponse, timingOffset] = channel(txWaveform, carrier);

    [noise, nPowerPerRE] = generateAWGNcomplex(SNRdB, nRxAnts, waveformInfo.Nfft, size(rxWaveform));
    rxWaveform = rxWaveform + noise;

    % ---- Timing synchronization ----
    if perfectEstimation
        offset = timingOffset;
    else
        [t,mag] = nrTimingEstimate(carrier, rxWaveform, dmrsIndices, dmrsSymbols);
        offset = skipWeakTimingOffset(offset, t, mag);
    end
    rxWaveform = rxWaveform(1+offset:end,:);

    % ---- OFDM demodulation ----
    rxGrid = nrOFDMDemodulate(carrier, rxWaveform);

    % ---- Channel estimation and precoder update ----
    if perfectEstimation
        estChGridAnts = ofdmChannelResponse;
        noiseEst = nPowerPerRE;

        % Compute next-slot precoder from antenna-domain estimate
        newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet, pdsch.NumLayers, estChGridAnts);

        % Convert antenna-domain estimate to layer-domain estimate for equalization
        estChGridLayers = precodeChannelEstimate(estChGridAnts, precodingWeights.');
    else
        % Practical estimate is layer->RX and includes precoding effect
        [estChGridLayers, noiseEst] = nrChannelEstimate(carrier, rxGrid, dmrsIndices, dmrsSymbols, ...
            "CDMLengths", pdsch.DMRS.CDMLengths);

        % Remove precoding effect to estimate antenna-domain channel
        estChGridAnts = precodeChannelEstimate(estChGridLayers, conj(precodingWeights));

        % Compute next-slot precoder
        newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet, pdsch.NumLayers, estChGridAnts);
    end

    % ---- Equalization ----
    [pdschRx, pdschHest] = nrExtractResources(pdschIndices, rxGrid, estChGridLayers);
    [pdschEq, csi] = nrEqualizeMMSE(pdschRx, pdschHest, noiseEst);

    % ---- Constellation plot ----
    constPlot.ChannelNames = "Layer " + (pdsch.NumLayers:-1:1);
    constPlot.ShowLegend = true;
    constPlot(fliplr(pdschEq));

    % ---- PDSCH decode to LLRs ----
    [dlschLLRs, rxSymbols] = nrPDSCHDecode(carrier, pdsch, pdschEq, noiseEst);

    % CSI scaling of LLRs
    csi = nrLayerDemap(csi);
    for cw = 1:pdsch.NumCodewords
        Qm = length(dlschLLRs{cw}) / length(rxSymbols{cw}); % bits per symbol
        csi{cw} = repmat(csi{cw}.', Qm, 1);
        dlschLLRs{cw} = dlschLLRs{cw} .* csi{cw}(:);
    end

    % ---- DL-SCH decode (uses soft buffers per process ID) ----
    decodeDLSCH.TransportBlockLength = trBlkSizes;
    [~, blkerr] = decodeDLSCH(dlschLLRs, pdsch.Modulation, pdsch.NumLayers, rv, pid0);

    % ---- HARQ state update ----
    harq = updateHARQManager(harq, pidx, blkerr);

    slotAnyBlkErr(nSlot+1) = any(blkerr);

    fprintf("Slot %2d | PID %2d | RV %s | BLKERR %s | NewData %s | TxN %s\n", ...
        nSlot, pid0, mat2str(rv), mat2str(blkerr), ...
        mat2str(harq.NewData(pidx,:)), mat2str(harq.TxNumber(pidx,:)));

end

fprintf("\nOverall BLER (any CW error) = %.4f over %d slots\n", mean(slotAnyBlkErr), totalNoSlots);

%% ====================== Local functions ======================

function harq = initHARQManager(Nproc, rvSeq, Ncw)
% Minimal HARQ manager replacing HARQEntity.
% State is stored in arrays sized [Nproc x Ncw] where processes are indexed 1..Nproc.
    harq.Nproc = Nproc;
    harq.Ncw   = Ncw;
    harq.rvSeq = rvSeq(:).';         % 1xNrv
    harq.Nrv   = numel(harq.rvSeq);

    harq.RVIdx = ones(Nproc,Ncw);    % index into rvSeq (starts at 1)
    harq.RedundancyVersion = repmat(harq.rvSeq(1), Nproc, Ncw);

    harq.TxNumber = ones(Nproc,Ncw); % # transmissions for current TB
    harq.NewData  = true(Nproc,Ncw); % need new TB initially
    harq.SequenceTimeout = false(Nproc,Ncw);
end

function harq = updateHARQManager(harq, pidx, blkerr)
% Update HARQ state for a given 1-based process index pidx after decoding.
    blkerr = blkerr(:).'; % row

    for cw = 1:harq.Ncw
        if ~blkerr(cw)
            % ACK: next time start new TB
            harq.NewData(pidx,cw) = true;
            harq.SequenceTimeout(pidx,cw) = false;

            harq.RVIdx(pidx,cw) = 1;
            harq.RedundancyVersion(pidx,cw) = harq.rvSeq(1);
            harq.TxNumber(pidx,cw) = 1;
        else
            % NACK: move to next RV if available else timeout -> new TB
            if harq.RVIdx(pidx,cw) < harq.Nrv
                harq.NewData(pidx,cw) = false; % retransmit same buffered TB
                harq.RVIdx(pidx,cw) = harq.RVIdx(pidx,cw) + 1;
                harq.RedundancyVersion(pidx,cw) = harq.rvSeq(harq.RVIdx(pidx,cw));
                harq.TxNumber(pidx,cw) = harq.TxNumber(pidx,cw) + 1;
                harq.SequenceTimeout(pidx,cw) = false;
            else
                harq.NewData(pidx,cw) = true;
                harq.SequenceTimeout(pidx,cw) = true;

                harq.RVIdx(pidx,cw) = 1;
                harq.RedundancyVersion(pidx,cw) = harq.rvSeq(1);
                harq.TxNumber(pidx,cw) = 1;
            end
        end
    end
end

function [noise,nPowerPerRE] = generateAWGNcomplex(SNRdB,nRxAnts,Nfft,sizeRxWaveform)
% Complex AWGN for SNR per RE per RX antenna.
    SNR = 10^(SNRdB/10);
    N0 = 1/sqrt(nRxAnts*double(Nfft)*SNR);
    w = (randn(sizeRxWaveform) + 1i*randn(sizeRxWaveform)) / sqrt(2);
    noise = N0 * w;
    nPowerPerRE = (N0^2) * double(Nfft);
end

function offsetOut = skipWeakTimingOffset(offsetPrev,t,mag)
% Minimal replacement for hSkipWeakTimingOffset.
% Keep previous offset if correlation is weak.

    if isempty(mag) || all(mag(:) == 0)
        offsetOut = offsetPrev;
        return;
    end

    peak   = max(mag(:));
    avgMag = mean(mag(:)) + eps;

    % Heuristic threshold: accept new timing only if peak is strong
    if peak/avgMag < 3
        offsetOut = offsetPrev;
    else
        offsetOut = t;
    end
end

function wtx = getPrecodingMatrix(PRBSet,NLayers,hestGrid)
% SVD-based precoder from antenna-domain channel estimate averaged over allocation.
    allocSc = (1:12)' + 12*PRBSet(:).';
    allocSc = allocSc(:);

    [~,~,R,P] = size(hestGrid);
    estAllocGrid = hestGrid(allocSc,:,:,:);
    Hest = permute(mean(reshape(estAllocGrid,[],R,P),1),[2 3 1]); % R x P

    [~,~,V] = svd(Hest);
    wtx = V(:,1:NLayers).';
    wtx = wtx/sqrt(NLayers);
end

function estChannelGrid = getInitialChannelEstimate(channel,carrier)
% Initial perfect antenna-domain channel estimate for precoder bootstrap.
    chClone = channel.clone();
    chClone.release();
    chClone.ChannelFiltering = false;
    chClone.ChannelResponseOutput = "ofdm-response";
    estChannelGrid = chClone(carrier);
end

function refPoints = getConstellationRefPoints(mod)
% Reference constellation points for comm.ConstellationDiagram
    switch mod
        case "QPSK"
            nPts = 4;
        case "16QAM"
            nPts = 16;
        case "64QAM"
            nPts = 64;
        case "256QAM"
            nPts = 256;
        otherwise
            error("Unsupported modulation: %s", string(mod));
    end
    binaryValues = int2bit(0:nPts-1,log2(nPts));
    refPoints = nrSymbolModulate(binaryValues(:),mod);
end

function estChannelGrid = precodeChannelEstimate(estChannelGrid,W)
% Apply precoding matrix W to the last dimension of the channel estimate.
    K = size(estChannelGrid,1);
    L = size(estChannelGrid,2);
    R = size(estChannelGrid,3);
    estChannelGrid = reshape(estChannelGrid,K*L*R,[]);
    estChannelGrid = estChannelGrid*W;
    estChannelGrid = reshape(estChannelGrid,K,L,R,[]);
end
