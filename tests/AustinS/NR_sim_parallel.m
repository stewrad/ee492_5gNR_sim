function result = NR_sim_parallel(cfg)
% NR_sim_parallel(cfg)
% cfg fields:
%   runIdx, SNRdB, Modulation, NHARQProcesses, rvSeq,
%   nTxAnts, nRxAnts, NumLayers, DelayProfile, logging, dispPlots, seed

%% ------------------------------------------------------------------------
% Unpack config
% -------------------------------------------------------------------------
SNRdB          = cfg.SNRdB;
Modulation     = cfg.Modulation;
NHARQProcesses = cfg.NHARQProcesses;
rvSeq          = cfg.rvSeq;
nTxAnts        = cfg.nTxAnts;
nRxAnts        = cfg.nRxAnts;
NumLayers      = cfg.NumLayers;
DelayProfile   = cfg.DelayProfile;
dispPlots      = cfg.dispPlots;
logging        = cfg.logging;

% Honor per-run seed from parfor sweep
rng(cfg.seed, 'twister');

%% ------------------------------------------------------------------------
% Result struct + safe defaults
% -------------------------------------------------------------------------
result = struct();
result.runIdx = cfg.runIdx;
result.logFile = "";
result.SNRdB = SNRdB;
result.Modulation = string(Modulation);
result.NHARQProcesses = NHARQProcesses;
result.DelayProfile = string(DelayProfile);

fid = -1;
logFile = "";
logBuffer = strings(0,1);
figFileConst = "";

%% ------------------------------------------------------------------------
% Logging / output setup
% -------------------------------------------------------------------------
if logging == 1
    logDir = fullfile(pwd, 'logs');
    if ~exist(logDir,'dir')
        mkdir(logDir);
    end

    figDir = fullfile(pwd,'figures');
    if ~exist(figDir,'dir')
        mkdir(figDir);
    end

    rvStr  = sprintf('%d', rvSeq);
    antStr = sprintf('%dx%d', nTxAnts, nRxAnts);

    ts = char(datetime('now','Format','yyyyMMdd_HHmmss'));

    logFile = fullfile(logDir, ...
        sprintf('%s_SNR%.1f_%s_NHARQ%d_rvSeq[%s]_%s_%s_run%05d.txt', ...
        Modulation, ...
        SNRdB, ...
        DelayProfile, ...
        NHARQProcesses, ...
        rvStr, ...
        antStr, ...
        ts, ...
        cfg.runIdx));

    fid = fopen(logFile, 'w');
    if fid < 0
        error('Could not open log file: %s', logFile);
    end

    [~, baseName, ~] = fileparts(logFile);
    figFileConst = fullfile(figDir, baseName + "_Constellation.png");

    result.logFile = string(logFile);
end

try
    %% --------------------------------------------------------------------
    % Specify SNR, number of slots, and perfect channel estimation flag
    % ---------------------------------------------------------------------
    totalNoSlots = 500;
    perfectEstimation = false;

    %% --------------------------------------------------------------------
    % Carrier configuration
    % ---------------------------------------------------------------------
    carrier = nrCarrierConfig;
    % note that nrCarrierConfig defaults to SCS = 15 kHz unless
    % carrier.SubcarrierSpacing is set explicitly elsehwere, e.g. 30 kHz
    carrier.SubcarrierSpacing = 15;   % kHz

    %% --------------------------------------------------------------------
    % PDSCH and DM-RS config
    % ---------------------------------------------------------------------
    pdsch = nrPDSCHConfig;
    pdsch.Modulation = Modulation;
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
        logBuffer(end+1) = sprintf('%15s: "%s"', 'Modulation', Modulation);
        logBuffer(end+1) = sprintf('%15s: %d', 'NHARQProcesses', NHARQProcesses);
        logBuffer(end+1) = sprintf('%15s: [%s]', 'rvSeq', strjoin(string(rvSeq), ' '));
        logBuffer(end+1) = sprintf('%15s: %d', 'nTxAnts', nTxAnts);
        logBuffer(end+1) = sprintf('%15s: %d', 'nRxAnts', nRxAnts);
        logBuffer(end+1) = sprintf('%15s: %d', 'NumLayers', NumLayers);
        logBuffer(end+1) = sprintf('%15s: "%s"', 'DelayProfile', DelayProfile);
        logBuffer(end+1) = "===========================";
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
    constPlot.ReferenceConstellation = getConstellationRefPoints(pdsch.Modulation);
    constPlot.EnableMeasurements = 1;

    offset = 0;

    totalTxBits = 0;
    totalRxBits = 0;
    totalTransmissions = 0;
    totalInitialTransmissions = 0;
    totalRetransmissions = 0;
    totalCodedBits = 0;
    blockErrors = 0;
    successfulBlocks = 0;
    transmissionAttempts = zeros(totalNoSlots, 1); %#ok<NASGU>

    perSlotSuccess = zeros(totalNoSlots, 1); %#ok<NASGU>
    perSlotBER = zeros(totalNoSlots, 1);

    estChannelGrid = getInitialChannelEstimate(channel,carrier);
    newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChannelGrid);

    instantaneousSINR = zeros(totalNoSlots, pdsch.NumLayers);
    avgChannelGain = zeros(totalNoSlots, 1);
    conditionNumber = zeros(totalNoSlots, 1);

    % EVM/MER arrays
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

    %% --------------------------------------------------------------------
    % Main simulation loop
    % ---------------------------------------------------------------------
    for nSlot = 0:totalNoSlots-1
        % attempt at correcting counting
        isNewDataCW = false(pdsch.NumCodewords,1);

        carrier.NSlot = nSlot;

        [pdschIndices,pdschInfo] = nrPDSCHIndices(carrier,pdsch);

        Xoh_PDSCH = 0;
        trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet),pdschInfo.NREPerPRB,codeRate,Xoh_PDSCH);

        for cwIdx = 1:pdsch.NumCodewords
        
            isNewDataCW(cwIdx) = harqEntity.NewData(cwIdx);

            if harqEntity.NewData(cwIdx)
                trBlk = randi([0 1],trBlkSizes(cwIdx),1);
                setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);

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
            H_matrix = squeeze(H_layer(size(H_layer,1)/2, :, :));
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
            M = evmMerMetricsDD(layerSyms, refConst_col);

            EVMrms_pct(nSlot+1, layerIdx) = M.RmsEVM_pct;
            EVMpk_pct (nSlot+1, layerIdx) = M.PeakEVM_pct;
            EVMrms_dB (nSlot+1, layerIdx) = M.AvgEVM_dB;
            EVMpk_dB  (nSlot+1, layerIdx) = M.PeakEVM_dB;
            MERavg_dB (nSlot+1, layerIdx) = M.AvgMER_dB;

            if logging == 1
                logBuffer(end+1) = sprintf(['Slot %3d | Layer %d | RMS EVM = %6.2f%% | Peak EVM = %6.2f%% | ' ...
                     'Avg EVM = %6.2f dB | Peak EVM = %6.2f dB | Avg MER = %6.2f dB'], ...
                     nSlot, layerIdx, M.RmsEVM_pct, M.PeakEVM_pct, M.AvgEVM_dB, M.PeakEVM_dB, M.AvgMER_dB);
            end
        end

        for layerIdx = 1:pdsch.NumLayers
            layerSymbols = pdschEq(:, layerIdx);
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
        [decbits,blkerr] = decodeDLSCH(dlschLLRs,pdsch.Modulation,pdsch.NumLayers, ...
            harqEntity.RedundancyVersion,harqEntity.HARQProcessID); %#ok<NASGU>

        % ATTEMPT correcting Throughput Efficiency calc
        for cwIdx = 1:pdsch.NumCodewords
            % Count every HARQ transmission attempt
            totalTransmissions = totalTransmissions + 1;
            % Count coded bits actually transmitted
            totalCodedBits = totalCodedBits + pdschInfo.G;
            % Count transmitted information bits on every attempt
            % (new transmission + retransmissions)
            totalTxBits = totalTxBits + trBlkSizes(cwIdx);

            if isNewDataCW(cwIdx)
                totalInitialTransmissions = totalInitialTransmissions + 1;
            else
                totalRetransmissions = totalRetransmissions + 1;
            end

            % Count successfully delivered information bits once per
            % successful decode event
            if ~blkerr(cwIdx)
                successfulBlocks = successfulBlocks + 1;
                totalRxBits = totalRxBits + trBlkSizes(cwIdx);
                perSlotSuccess(nSlot+1) = 1;
            end
        end

        statusReport = updateAndAdvance(harqEntity,blkerr,trBlkSizes,pdschInfo.G);

        [perfState, perfLog(nSlot+1)] = logPerfMetrics(perfState, nSlot, trBlkSizes, blkerr, harqEntity, slotDuration_s);

        if logging == 1
            logBuffer(end+1) = sprintf(['Slot %3d | NewTB=%d | TxBits=%6d | RxBits=%6d | ' ...
                 'AttemptBLER=%.3f | CumAttemptBLER=%.3f | ' ...
                 'FinalBLER=%.3f | Eff=%.1f%% | InstThr=%.2f | AvgThr=%.2f'], ...
                nSlot, perfLog(nSlot+1).isNewTB, perfLog(nSlot+1).txBitsThisSlot, ...
                perfLog(nSlot+1).rxBitsThisSlot, ...
                perfLog(nSlot+1).attemptBLER_thisSlot, perfLog(nSlot+1).attemptBLER_cum, ...
                perfLog(nSlot+1).finalBLER, perfLog(nSlot+1).throughputEfficiencyPct, ...
                perfLog(nSlot+1).instThroughputMbps, perfLog(nSlot+1).avgThroughputMbps);
        end

        % totalTransmissions = totalTransmissions + 1;
        % 
        % for cwIdx = 1:pdsch.NumCodewords
        %     if harqEntity.NewData(cwIdx)
        %         totalInitialTransmissions = totalInitialTransmissions + 1;
        %         totalTxBits = totalTxBits + trBlkSizes(cwIdx);
        %     else
        %         totalRetransmissions = totalRetransmissions + 1;
        %     end
        % 
        %     if ~blkerr(cwIdx)
        %         successfulBlocks = successfulBlocks + 1;
        %         totalRxBits = totalRxBits + trBlkSizes(cwIdx);
        %         perSlotSuccess(nSlot+1) = 1;
        %     else
        %         if harqEntity.SequenceTimeout(cwIdx)
        %             blockErrors = blockErrors + 1;
        %         end
        %     end
        % end

        perSlotBER(nSlot+1) = sum(blkerr) / pdsch.NumCodewords;

        if logging == 1
            logBuffer(end+1) = sprintf('Slot %d. %s', nSlot, statusReport);
        end

        lastEstChGridLayers = estChGridLayers; %#ok<NASGU>
        lastPdschEq         = pdschEq;
    end

    %% --------------------------------------------------------------------
    % Save constellation figure
    % ---------------------------------------------------------------------
    refConst   = getConstellationRefPoints(pdsch.Modulation);
    refConst   = refConst(:);
    layerColors = [0.15 0.45 0.85;
                   0.85 0.33 0.10;
                   0.47 0.67 0.19;
                   0.64 0.08 0.18];

    figConst = figure('Name','Constellation','Visible','off');
    figConst.Position = [100 100 580 560];
    set(figConst, 'Color', 'white');
    ax = axes('Parent', figConst);
    set(ax, 'Color', 'white', 'XColor', 'black', 'YColor', 'black', ...
            'GridColor', [0.82 0.82 0.82], 'FontSize', 9);
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
    xlabel(ax, 'In-phase Amplitude', 'Color', 'black');
    ylabel(ax, 'Quadrature Amplitude', 'Color', 'black');
    title(ax, sprintf('Equalized PDSCH Constellation  |  %s  |  SNR = %.1f dB  |  Last Slot', ...
          pdsch.Modulation, SNRdB), 'Color', 'black');

    lgd = legend(ax, 'Location', 'northeast', 'FontSize', 8);
    set(lgd, 'Color', 'white', 'EdgeColor', 'black', 'TextColor', 'black');

    ann_lines = {};
    for lyr = 1:pdsch.NumLayers
        rmsEVM  = EVMrms_pct(end, lyr);
        pkEVM   = EVMpk_pct (end, lyr);
        evmDB   = EVMrms_dB (end, lyr);
        pkEvmDB = EVMpk_dB  (end, lyr);
        merAvg  = MERavg_dB (end, lyr);
        ann_lines{end+1} = sprintf('--- Layer %d ---', lyr); %#ok<SAGROW>
        ann_lines{end+1} = sprintf('RMS EVM  : %6.2f %%', rmsEVM); %#ok<SAGROW>
        ann_lines{end+1} = sprintf('Peak EVM : %6.2f %%', pkEVM); %#ok<SAGROW>
        ann_lines{end+1} = sprintf('Avg EVM  : %6.2f dB', evmDB); %#ok<SAGROW>
        ann_lines{end+1} = sprintf('Peak EVM : %6.2f dB', pkEvmDB); %#ok<SAGROW>
        ann_lines{end+1} = sprintf('Avg MER  : %6.2f dB', merAvg); %#ok<SAGROW>
        if lyr < pdsch.NumLayers
            ann_lines{end+1} = ''; %#ok<SAGROW>
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
        logBuffer(end+1) = sprintf('Constellation figure saved to: %s\n', figFileConst);
    else
        exportgraphics(figConst, fullfile(pwd, 'Constellation.png'), 'Resolution', 150);
    end
    close(figConst);

    %% --------------------------------------------------------------------
    % Summary logging
    % ---------------------------------------------------------------------
    logBuffer(end+1) = sprintf('========================================');
    logBuffer(end+1) = sprintf('     5G NR HARQ Simulation Results      ');
    logBuffer(end+1) = sprintf('========================================');

    logBuffer(end+1) = sprintf('\n--- Configuration ---');
    logBuffer(end+1) = sprintf('HARQ Type: Chase Combining');
    logBuffer(end+1) = sprintf('RV Sequence: [%s]', strjoin(string(rvSeq), ' '));
    logBuffer(end+1) = sprintf('Number of HARQ Processes: %d', NHARQProcesses);
    logBuffer(end+1) = sprintf('Modulation: %s', pdsch.Modulation);
    logBuffer(end+1) = sprintf('Number of Layers: %d', pdsch.NumLayers);
    logBuffer(end+1) = sprintf('Code Rate: %.4f', codeRate);
    logBuffer(end+1) = sprintf('SNR: %.2f dB', SNRdB);
    logBuffer(end+1) = sprintf('Number of Slots: %d', totalNoSlots);

    logBuffer(end+1) = sprintf('\n--- Transmission Statistics ---');
    logBuffer(end+1) = sprintf('Total Transmissions: %d', totalTransmissions);
    logBuffer(end+1) = sprintf('Initial Transmissions: %d', totalInitialTransmissions);
    logBuffer(end+1) = sprintf('Retransmissions: %d', totalRetransmissions);
    logBuffer(end+1) = sprintf('Average Transmissions per TB: %.2f', ...
        totalTransmissions / max(1,totalInitialTransmissions));
    logBuffer(end+1) = sprintf('Retransmission Rate: %.2f%%', ...
        (totalRetransmissions / max(1,totalInitialTransmissions + totalRetransmissions)) * 100);

    finalBLER = (totalInitialTransmissions - successfulBlocks) / max(1,totalInitialTransmissions);
    attemptBLER = perfState.attemptsFailed / max(1, perfState.attemptsTotal);

    logBuffer(end+1) = sprintf('\n--- Performance Metrics ---');
    logBuffer(end+1) = sprintf('Successful Blocks: %d / %d', successfulBlocks, totalInitialTransmissions);
    logBuffer(end+1) = sprintf('Failed Blocks (after max retx): %d', blockErrors);
    logBuffer(end+1) = sprintf('Attempt BLER (pre-HARQ): %.4f (%.2f%%)', attemptBLER, attemptBLER*100);
    logBuffer(end+1) = sprintf('Final BLER (post-HARQ):  %.4f (%.2f%%)', finalBLER, finalBLER*100);

    logBuffer(end+1) = sprintf('\n--- Throughput Analysis ---');
    % logBuffer(end+1) = sprintf('Total Bits Transmitted: %d', totalTxBits);
    % logBuffer(end+1) = sprintf('Total Bits Received (success): %d', totalRxBits);
    % logBuffer(end+1) = sprintf('Throughput Efficiency: %.2f%%', (totalRxBits / max(1,totalTxBits)) * 100);
    % logBuffer(end+1) = sprintf('Effective Code Rate: %.4f', totalRxBits / max(1,(totalTransmissions * pdschInfo.G)));
    logBuffer(end+1) = sprintf('Total Information Bits Sent (all HARQ attempts): %d', totalTxBits);
    logBuffer(end+1) = sprintf('Total Information Bits Delivered Successfully: %d', totalRxBits);
    logBuffer(end+1) = sprintf('Throughput Efficiency: %.2f%%', ...
        (totalRxBits / max(1,totalTxBits)) * 100);
    logBuffer(end+1) = sprintf('Effective Throughput Efficiency (InfoBits/CodedBits): %.4f', ...
    totalRxBits / max(1,totalCodedBits));
    logBuffer(end+1) = sprintf('Nominal Code Rate (per TB): %.4f', ...
    trBlkSizes(1) / pdschInfo.G);

    avgBitsPerSlot = totalRxBits / totalNoSlots;
    slotDuration = carrier.SlotsPerSubframe / (carrier.SubcarrierSpacing / 15e3) * 1e-3;
    throughput_Mbps = (avgBitsPerSlot / slotDuration) / 1e6;

    logBuffer(end+1) = sprintf('Average Throughput: %.2f Mbps', throughput_Mbps);
    logBuffer(end+1) = sprintf('Average Bits per Slot: %.0f bits', avgBitsPerSlot);

    latency_est = (1 + (totalRetransmissions/totalInitialTransmissions)) * slotDuration_s;
    logBuffer(end+1) = sprintf('\n--- Latency Calculation ---');
    logBuffer(end+1) = sprintf('Time slot Duration (ms): %.6f', slotDuration_s * 1e3);
    logBuffer(end+1) = sprintf('Estimated Latency (ms): %.6f', latency_est * 1e3);

    logBuffer(end+1) = sprintf('\n--- Spectral Efficiency ---');
    numRBs = numel(pdsch.PRBSet);
    bandwidth_Hz = numRBs * 12 * carrier.SubcarrierSpacing * 1e3;
    spectralEfficiency = (avgBitsPerSlot / slotDuration) / bandwidth_Hz;
    logBuffer(end+1) = sprintf('Spectral Efficiency: %.4f bits/s/Hz', spectralEfficiency);

    logBuffer(end+1) = sprintf('\n========================================\n');
    logBuffer(end+1) = sprintf('--- Channel Model Details ---');
    logBuffer(end+1) = sprintf('Delay Profile: %s', channel.DelayProfile);
    logBuffer(end+1) = sprintf('RMS Delay Spread: %.2f ns', rmsDelaySpread * 1e9);
    logBuffer(end+1) = sprintf('Mean Excess Delay: %.2f ns', meanDelay * 1e9);
    logBuffer(end+1) = sprintf('Maximum Channel Delay: %d samples (%.2f μs)', ...
        chInfo.MaximumChannelDelay, ...
        chInfo.MaximumChannelDelay / channel.SampleRate * 1e6);

    logBuffer(end+1) = sprintf('\n--- Channel Quality Statistics ---');
    logBuffer(end+1) = sprintf('Average Channel Gain: %.2f dB', mean(avgChannelGain));
    logBuffer(end+1) = sprintf('Channel Gain Std Dev: %.2f dB', std(avgChannelGain));
    logBuffer(end+1) = sprintf('Average Condition Number: %.2f', mean(conditionNumber(conditionNumber > 0)));

    logBuffer(end+1) = sprintf('\n--- SINR Statistics (per layer) ---');
    for layerIdx = 1:pdsch.NumLayers
        logBuffer(end+1) = sprintf('Layer %d - Mean SINR: %.2f dB, Std: %.2f dB', ...
            layerIdx, ...
            mean(instantaneousSINR(:, layerIdx)), ...
            std(instantaneousSINR(:, layerIdx)));
    end

    coherenceBW_Hz = 1 / (5 * rmsDelaySpread);
    logBuffer(end+1) = sprintf('\nCoherence Bandwidth (50%% corr): %.2f kHz', coherenceBW_Hz / 1e3);

    if isprop(channel, 'MaximumDopplerShift') && channel.MaximumDopplerShift > 0
        coherenceTime_s = 9 / (16 * pi * channel.MaximumDopplerShift);
        logBuffer(end+1) = sprintf('Coherence Time (50%% corr): %.2f ms', coherenceTime_s * 1e3);
    end

    %% --------------------------------------------------------------------
    % Final log flush
    % ---------------------------------------------------------------------
    if logging == 1 && fid >= 0
        fprintf(fid, '%s\n', logBuffer);
        fprintf(fid, '\nRun End: %s\n', datestr(now));
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