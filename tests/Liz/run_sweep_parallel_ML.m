% run_sweep_parallel_ML.m
% Parallel parameter sweep for NR_sim_parallel(cfg)
% using ML-selected modulation from an Excel decision table

clear all;
close all;
clc;

sweepTimer = tic;

%% ML decision table file
mlFile = 'ml_results_output_extended.xlsx';   % change if needed
mlSheet = 'Decision_Table';

% Read table and preserve original Excel headers
mlTable = readtable(mlFile, 'Sheet', mlSheet, 'VariableNamingRule', 'preserve');

% Check required columns exist
requiredCols = ["SNRdB", "NHARQProcesses", "Recommended Modulation"];
for k = 1:numel(requiredCols)
    if ~ismember(requiredCols(k), string(mlTable.Properties.VariableNames))
        error('Missing required column "%s" in %s (%s sheet).', ...
            requiredCols(k), mlFile, mlSheet);
    end
end

%% Sweep parameters
SNRdB_vals          = [0 2 4 6 8 10 12 14 16];
NHARQProcesses_vals = [1 4 8];
rvSeq_vals          = {[0 2 3 1]};
nTxAnts_vals        = [2];
nRxAnts_vals        = [2];
DelayProfile_vals   = {"TDL-C"};
NumLayers           = 2;

%% Build all parameter combinations
[gS,gH,gRV,gTx,gRx,gDP] = ndgrid( ...
    1:numel(SNRdB_vals), ...
    1:numel(NHARQProcesses_vals), ...
    1:numel(rvSeq_vals), ...
    1:numel(nTxAnts_vals), ...
    1:numel(nRxAnts_vals), ...
    1:numel(DelayProfile_vals));

idxList = [gS(:), gH(:), gRV(:), gTx(:), gRx(:), gDP(:)];
totalRuns = size(idxList, 1);

fprintf('\n========================================================\n');
fprintf('  Starting ML-based parallel sweep\n');
fprintf('  Total combinations: %d\n', totalRuns);
fprintf('  Decision table: %s (%s)\n', mlFile, mlSheet);
fprintf('  Start time: %s\n', datetime('now'));
fprintf('========================================================\n\n');

%% Start or reset parallel pool
numWorkers = 2;

poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local', numWorkers);
elseif poolobj.NumWorkers ~= numWorkers
    delete(poolobj);
    parpool('local', numWorkers);
end

%% Preallocate results
results(totalRuns,1) = struct( ...
    'runIdx', 0, ...
    'ok', false, ...
    'skipped', false, ...
    'errorMessage', "", ...
    'logFile', "", ...
    'SNRdB', 0, ...
    'Modulation', "", ...
    'NHARQProcesses', 0, ...
    'rvSeq', [], ...
    'nTxAnts', 0, ...
    'nRxAnts', 0, ...
    'NumLayers', 0, ...
    'DelayProfile', "" );

%% Progress tracking
dq = parallel.pool.DataQueue;
afterEach(dq, @updateProgress);

%% Parallel sweep
parfor runIdx = 1:totalRuns

    cfg = struct();
    cfg.runIdx          = runIdx;
    cfg.SNRdB           = SNRdB_vals(idxList(runIdx,1));
    cfg.NHARQProcesses  = NHARQProcesses_vals(idxList(runIdx,2));
    cfg.rvSeq           = rvSeq_vals{idxList(runIdx,3)};
    cfg.nTxAnts         = nTxAnts_vals(idxList(runIdx,4));
    cfg.nRxAnts         = nRxAnts_vals(idxList(runIdx,5));
    cfg.DelayProfile    = DelayProfile_vals{idxList(runIdx,6)};
    cfg.NumLayers       = NumLayers;
    cfg.logging         = 1;
    cfg.dispPlots       = 0;
    cfg.seed            = 200000 + runIdx;

    % Look up ML-selected modulation
    match = mlTable.("SNRdB") == cfg.SNRdB & ...
            mlTable.("NHARQProcesses") == cfg.NHARQProcesses;

    r = struct( ...
        'runIdx', runIdx, ...
        'ok', false, ...
        'skipped', false, ...
        'errorMessage', "", ...
        'logFile', "", ...
        'SNRdB', cfg.SNRdB, ...
        'Modulation', "", ...
        'NHARQProcesses', cfg.NHARQProcesses, ...
        'rvSeq', cfg.rvSeq, ...
        'nTxAnts', cfg.nTxAnts, ...
        'nRxAnts', cfg.nRxAnts, ...
        'NumLayers', cfg.NumLayers, ...
        'DelayProfile', string(cfg.DelayProfile) );

    try
        if sum(match) ~= 1
            error('Decision table lookup failed for SNR=%g, NHARQ=%d.', ...
                cfg.SNRdB, cfg.NHARQProcesses);
        end

        cfg.Modulation = char(string(mlTable.("Recommended Modulation")(match)));
        r.Modulation = string(cfg.Modulation);

        % Skip invalid antenna/layer combinations if needed
        if (cfg.nTxAnts < cfg.NumLayers) || (cfg.nRxAnts < cfg.NumLayers)
            r.skipped = true;
            r.ok = false;
            r.errorMessage = "Invalid antenna/layer combo";
        else
            out = NR_sim_parallel(cfg);
            r.ok = true;

            if isfield(out, 'logFile')
                r.logFile = string(out.logFile);
            end
        end

    catch ME
        r.ok = false;
        r.skipped = false;
        r.errorMessage = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
    end

    results(runIdx) = r;
    send(dq, 1);
end

%% Summary
finished = sum([results.ok]);
skipped  = sum([results.skipped]);
failed   = totalRuns - finished - skipped;

save('sweep_results_parallel_ML.mat', 'results');

% Print first failure if any
badIdx = find(~[results.ok] & ~[results.skipped], 1, 'first');
if ~isempty(badIdx)
    fprintf('\nFirst runtime failure was run %d\n', badIdx);
    fprintf('%s\n', results(badIdx).errorMessage);
end

fprintf('\n========================================================\n');
fprintf('  ML sweep complete.\n');
fprintf('  Finished : %d\n', finished);
fprintf('  Skipped  : %d  (invalid antenna/layer combos)\n', skipped);
fprintf('  Failed   : %d  (runtime errors)\n', failed);
fprintf('  Total time: %s\n', formatDuration(toc(sweepTimer)));
fprintf('========================================================\n');

%% Local functions
function updateProgress(~)
    persistent count lastTime startTime

    if isempty(count)
        count = 0;
        startTime = tic;
        lastTime  = tic;
    end

    count = count + 1;

    if mod(count,10) == 0
        chunkTime = toc(lastTime);
        totalTime = toc(startTime);

        fprintf('Completed %d runs... | Last 10: %.2fs | Total: %s\n', ...
            count, chunkTime, formatDuration(totalTime));

        lastTime = tic;
    end
end

function str = formatDuration(seconds)
    h = floor(seconds/3600);
    m = floor(mod(seconds,3600)/60);
    s = mod(seconds,60);
    str = sprintf('%02dh:%02dm:%05.2fs', h, m, s);
end