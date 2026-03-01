% run_sweep.m
% Automated parameter sweep runner for nr_dlsch_tx_rx_updated.m
% Builds every combination of the values below and executes the main script
% for each one. Logs are written per-run (controlled by logging flag inside
% the main script). Progress is printed to the Command Window.
%
% WARNING: with the default value sets the full grid is 26,880 combinations.
% Narrow the sweep arrays before running to keep runtime manageable —
% e.g. pick a single DelayProfile and 2-3 SNR points to start.
%
% USAGE: run this file directly — do NOT run nr_dlsch_tx_rx_updated.m manually.
%   >> run_sweep
%
% NOTE: the main script must be on the MATLAB path (or in the same folder).
% -------------------------------------------------------------------------

% =========================================================================
% DEFINE SWEEP SPACE — edit these cell/numeric arrays to taste
% =========================================================================
SNRdB_vals        = [0 2 4 6 8 10 12];
Modulation_vals   = {"QPSK", "16QAM", "64QAM", "256QAM", "1024QAM"};   
NHARQProcesses_vals = [1 2 4 8 16];
rvSeq_vals        = {[0 2 3 1], [0 1 2 3], [0 3 2 1], [0 2 1 3]};
nTxAnts_vals      = [2 8];
nRxAnts_vals      = [2 8];
DelayProfile_vals = {"TDL-A", "TDL-B", "TDL-C"};

% NumLayers is derived automatically from antenna counts (see rule below).
% Override by setting a fixed value instead, e.g. NumLayers_fixed = 2;
% NumLayers_fixed   = [];   % leave empty [] to auto-derive
NumLayers = 2;

% =========================================================================
% BUILD COMBINATION LIST
% =========================================================================
% Use ndgrid to enumerate all index combinations, then linearise.
nS  = numel(SNRdB_vals);
nM  = numel(Modulation_vals);
nH  = numel(NHARQProcesses_vals);
nRV = numel(rvSeq_vals);
nTx = numel(nTxAnts_vals);
nRx = numel(nRxAnts_vals);
nDP = numel(DelayProfile_vals);

[gS, gM, gH, gRV, gTx, gRx, gDP] = ndgrid( ...
    1:nS, 1:nM, 1:nH, 1:nRV, 1:nTx, 1:nRx, 1:nDP);

idxList = [gS(:), gM(:), gH(:), gRV(:), gTx(:), gRx(:), gDP(:)];
totalRuns = size(idxList, 1);

fprintf('\n========================================================\n');
fprintf('  Parameter sweep: %d total combinations\n', totalRuns);
fprintf('========================================================\n\n');

% =========================================================================
% SWEEP LOOP
% =========================================================================
skipped  = 0;
finished = 0;
failed   = 0;
sweepTimer = tic;   % overall sweep wall-clock timer

for runIdx = 1:totalRuns

    % --- Unpack indices into parameter values ---
    SNRdB          = SNRdB_vals(        idxList(runIdx,1));
    Modulation     = Modulation_vals{   idxList(runIdx,2)};
    NHARQProcesses = NHARQProcesses_vals(idxList(runIdx,3));
    rvSeq          = rvSeq_vals{        idxList(runIdx,4)};
    nTxAnts        = nTxAnts_vals(      idxList(runIdx,5));
    nRxAnts        = nRxAnts_vals(      idxList(runIdx,6));
    DelayProfile   = DelayProfile_vals{ idxList(runIdx,7)};

    % --- Print progress ---
    rvStr = sprintf('%d', rvSeq);
    fprintf('[%4d/%d] SNR=%2gdB | %s | HARQ=%2d | rv=[%s] | %dx%d | L=%d | %s\n', ...
        runIdx, totalRuns, SNRdB, Modulation, NHARQProcesses, ...
        rvStr, nTxAnts, nRxAnts, NumLayers, DelayProfile);

    % --- Run the main script (all workspace vars are visible to it) ---
    runTimer = tic;
    try
        run('nr_dlsch_tx_rx_TEST.m');
        elapsed = toc(runTimer);
        finished = finished + 1;

        fprintf('\n\n\n         Done in %.1f s  |  elapsed total: %s  |  est. remaining: %s\n\n\n\n', ...
            elapsed, ...
            formatDuration(toc(sweepTimer)), ...
            formatDuration((toc(sweepTimer) / runIdx) * (totalRuns - runIdx)));

    catch ME
        elapsed = toc(runTimer);
        fprintf('  *** ERROR in run %d (%.1f s): %s ***\n', runIdx, elapsed, ME.message);
        failed = failed + 1;
    end

end % for runIdx

% =========================================================================
% SUMMARY
% =========================================================================
fprintf('\n========================================================\n');
fprintf('  Sweep complete.\n');
fprintf('  Finished : %d\n', finished);
fprintf('  Skipped  : %d  (invalid antenna/layer combos)\n', skipped);
fprintf('  Failed   : %d  (runtime errors)\n', failed);
fprintf('  Total time: %s\n', formatDuration(toc(sweepTimer)));
fprintf('========================================================\n');



function str = formatDuration(seconds)
% Convert a duration in seconds to a human-readable string: e.g. "1h 23m 45s"
    h = floor(seconds / 3600);
    m = floor(mod(seconds, 3600) / 60);
    s = floor(mod(seconds, 60));
    if h > 0
        str = sprintf('%dh %02dm %02ds', h, m, s);
    elseif m > 0
        str = sprintf('%dm %02ds', m, s);
    else
        str = sprintf('%ds', s);
    end
end