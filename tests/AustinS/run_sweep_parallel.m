% run_sweep_parallel.m
% Parallel parameter sweep for NR_sim_parallel(cfg)

clear all;
close all;
clc;

sweepTimer = tic;

% %% Sweep parameters
% SNRdB_vals          = [0 2 4 6 8 10 12];
% Modulation_vals     = {"QPSK","16QAM","64QAM","256QAM","1024QAM"};
% NHARQProcesses_vals = [1];  %[1 2 4 8 16];
% rvSeq_vals          = {[0 2 3 1], [0 0 0 0]}; %, [0 1 2 3], [0 3 2 1], [0 2 1 3], [0 0 0 0]};
% nTxAnts_vals        = [2 8];
% nRxAnts_vals        = [2 8];
% DelayProfile_vals   = {"TDL-C"}; %{"TDL-A","TDL-B","TDL-C"};
% NumLayers           = 2;

% LATENCY RUNS: 
SNRdB_vals          = [0 2 4 6 8 10 12 14 16 18];
Modulation_vals     = {"QPSK","16QAM","64QAM","256QAM","1024QAM"}; %{"256QAM"};
NHARQProcesses_vals = [1 4 16]; %[8 12 16 20 24]; 
rvSeq_vals          = {[0 2 3 1]}; %, [0 1 2 3], [0 3 2 1], [0 2 1 3], [0 0 0 0]};
nTxAnts_vals        = [2];
nRxAnts_vals        = [2];
DelayProfile_vals   = {"TDL-C"}; %{"TDL-A","TDL-B","TDL-C"};
NumLayers           = 2;


%% Build all parameter combinations
[gS,gM,gH,gRV,gTx,gRx,gDP] = ndgrid( ...
    1:numel(SNRdB_vals), ...
    1:numel(Modulation_vals), ...
    1:numel(NHARQProcesses_vals), ...
    1:numel(rvSeq_vals), ...
    1:numel(nTxAnts_vals), ...
    1:numel(nRxAnts_vals), ...
    1:numel(DelayProfile_vals));

idxList = [gS(:), gM(:), gH(:), gRV(:), gTx(:), gRx(:), gDP(:)];
totalRuns = size(idxList, 1);

fprintf('\n========================================================\n');
fprintf('  Starting parallel sweep\n');
fprintf('  Total combinations: %d\n', totalRuns);
fprintf('  Start time: %s\n', datetime('now'));
fprintf('========================================================\n\n');

%% Start or reset parallel pool
numWorkers = 10;   % change to 3 if desired

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
progressCounter = 0;
afterEach(dq, @updateProgress);

%% Parallel sweep
parfor runIdx = 1:totalRuns

    cfg = struct();
    cfg.runIdx          = runIdx;
    cfg.SNRdB           = SNRdB_vals(idxList(runIdx,1));
    cfg.Modulation      = Modulation_vals{idxList(runIdx,2)};
    cfg.NHARQProcesses  = NHARQProcesses_vals(idxList(runIdx,3));
    cfg.rvSeq           = rvSeq_vals{idxList(runIdx,4)};
    cfg.nTxAnts         = nTxAnts_vals(idxList(runIdx,5));
    cfg.nRxAnts         = nRxAnts_vals(idxList(runIdx,6));
    cfg.DelayProfile    = DelayProfile_vals{idxList(runIdx,7)};
    cfg.NumLayers       = NumLayers;
    cfg.logging         = 1;
    cfg.dispPlots       = 0;
    cfg.seed            = 100000 + runIdx;

    r = struct( ...
        'runIdx', runIdx, ...
        'ok', false, ...
        'skipped', false, ...
        'errorMessage', "", ...
        'logFile', "", ...
        'SNRdB', cfg.SNRdB, ...
        'Modulation', string(cfg.Modulation), ...
        'NHARQProcesses', cfg.NHARQProcesses, ...
        'rvSeq', cfg.rvSeq, ...
        'nTxAnts', cfg.nTxAnts, ...
        'nRxAnts', cfg.nRxAnts, ...
        'NumLayers', cfg.NumLayers, ...
        'DelayProfile', string(cfg.DelayProfile) );

    try
        % Skip invalid cases if needed
        if (cfg.nTxAnts < cfg.NumLayers) || (cfg.nRxAnts < cfg.NumLayers)
            r.skipped = true;
            r.ok = false;
            r.errorMessage = "Invalid antenna/layer combo";
        else
            out = NR_sim_parallel(cfg);
            % out = NR_sim_parallel_FBdelay(cfg);

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

save('sweep_results_parallel.mat', 'results');

% Print first failure if any
badIdx = find(~[results.ok] & ~[results.skipped], 1, 'first');
if ~isempty(badIdx)
    fprintf('\nFirst runtime failure was run %d\n', badIdx);
    fprintf('%s\n', results(badIdx).errorMessage);
end

fprintf('\n========================================================\n');
fprintf('  Sweep complete.\n');
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
        startTime = tic;      % total sweep timer
        lastTime  = tic;      % timer for 25-run chunks
    end

    count = count + 1;

    if mod(count,25) == 0
        chunkTime = toc(lastTime);
        totalTime = toc(startTime);

        fprintf('Completed %d runs... | Last 25: %.2fs | Total: %s\n', ...
            count, chunkTime, formatDuration(totalTime));

        lastTime = tic;   % reset chunk timer
    end
end

function str = formatDuration(seconds)
    h = floor(seconds/3600);
    m = floor(mod(seconds,3600)/60);
    s = mod(seconds,60);
    str = sprintf('%02dh:%02dm:%05.2fs', h, m, s);
end