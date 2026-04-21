
% run_sweep_parallel_ML.m
% Two-phase workflow:
%   1) Generate fixed-modulation dataset (16QAM and 64QAM)
%   2) Train a simple logistic-regression controller
%   3) Run ML-controlled sweep using the trained model
%
% Outputs:
%   - logs/...
%   - ml_dataset/... per-run dataset MAT files
%   - ml_model/ml_logistic_model.mat
%   - ml_model/ml_training_summary.txt
%   - sweep_results_parallel_ML.mat

clear all;
close all;
clc;

sweepTimer = tic;

%% ------------------------------------------------------------------------
% User settings
% -------------------------------------------------------------------------
SNRdB_vals          = 0:2:16;
NHARQProcesses_vals = [1 4 16];
rvSeq_vals          = {[0 2 3 1]};
nTxAnts_vals        = [2];
nRxAnts_vals        = [2];
DelayProfile_vals   = {"TDL-C"};
NumLayers           = 2;

fixedTrainingMods   = {"16QAM","64QAM"};
mlRiskThreshold     = 0.35;
numWorkers          = 10;
totalNoSlots        = 5000;

fprintf('\n========================================================\n');
fprintf('  Starting ML-assisted HARQ workflow\n');
fprintf('  Phase 1: dataset generation (fixed 16QAM / 64QAM)\n');
fprintf('  Phase 2: train logistic model\n');
fprintf('  Phase 3: ML-controlled sweep\n');
fprintf('  Start time: %s\n', datetime('now'));
fprintf('========================================================\n\n');

%% ------------------------------------------------------------------------
% Start/reset parallel pool
% -------------------------------------------------------------------------
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local', numWorkers);
elseif poolobj.NumWorkers ~= numWorkers
    delete(poolobj);
    parpool('local', numWorkers);
end

%% ------------------------------------------------------------------------
% Build base parameter combinations
% -------------------------------------------------------------------------
[gS,gH,gRV,gTx,gRx,gDP] = ndgrid( ...
    1:numel(SNRdB_vals), ...
    1:numel(NHARQProcesses_vals), ...
    1:numel(rvSeq_vals), ...
    1:numel(nTxAnts_vals), ...
    1:numel(nRxAnts_vals), ...
    1:numel(DelayProfile_vals));

baseIdxList = [gS(:), gH(:), gRV(:), gTx(:), gRx(:), gDP(:)];

%% ------------------------------------------------------------------------
% Phase 1: fixed-modulation dataset generation
% -------------------------------------------------------------------------
datasetCfgList = struct( ...
    'runIdx', {}, ...
    'SNRdB', {}, ...
    'NHARQProcesses', {}, ...
    'rvSeq', {}, ...
    'nTxAnts', {}, ...
    'nRxAnts', {}, ...
    'DelayProfile', {}, ...
    'NumLayers', {}, ...
    'logging', {}, ...
    'dispPlots', {}, ...
    'seed', {}, ...
    'totalNoSlots', {}, ...
    'controlMode', {}, ...
    'fixedModulation', {}, ...
    'datasetMode', {}, ...
    'mlRiskThreshold', {} );
runIdx = 0;

for m = 1:numel(fixedTrainingMods)
    for k = 1:size(baseIdxList,1)
        runIdx = runIdx + 1;
        cfg = struct();
        cfg.runIdx          = runIdx;
        cfg.SNRdB           = SNRdB_vals(baseIdxList(k,1));
        cfg.NHARQProcesses  = NHARQProcesses_vals(baseIdxList(k,2));
        cfg.rvSeq           = rvSeq_vals{baseIdxList(k,3)};
        cfg.nTxAnts         = nTxAnts_vals(baseIdxList(k,4));
        cfg.nRxAnts         = nRxAnts_vals(baseIdxList(k,5));
        cfg.DelayProfile    = DelayProfile_vals{baseIdxList(k,6)};
        cfg.NumLayers       = NumLayers;
        cfg.logging         = 1;
        cfg.dispPlots       = 0;
        cfg.seed            = 200000 + runIdx;
        cfg.totalNoSlots    = totalNoSlots;
        cfg.controlMode     = 'fixed';
        cfg.fixedModulation = fixedTrainingMods{m};
        cfg.datasetMode     = true;
        cfg.mlRiskThreshold = mlRiskThreshold;
        datasetCfgList(end+1) = cfg; %#ok<SAGROW>
    end
end

fprintf('Phase 1 dataset runs: %d\n', numel(datasetCfgList));
datasetResults = executeSweep(datasetCfgList, numWorkers, 'Phase 1 / Dataset');

datasetFinished = sum([datasetResults.ok]);
datasetFailed   = sum(~[datasetResults.ok] & ~[datasetResults.skipped]);
if datasetFinished == 0
    badIdx = find(~[datasetResults.ok] & ~[datasetResults.skipped], 1, 'first');
    if ~isempty(badIdx)
        error("Dataset generation failed before training. First failure:\n%s", datasetResults(badIdx).errorMessage);
    else
        error("Dataset generation produced no successful runs.");
    end
end
if datasetFailed > 0
    fprintf('Warning: %d dataset-generation runs failed. Training will use successful runs only.\n', datasetFailed);
end

%% ------------------------------------------------------------------------
% Combine dataset files and train logistic model
% -------------------------------------------------------------------------
fprintf('\nTraining logistic model from generated dataset...\n');
[mlModel, trainSummary] = trainMLModelFromResults(datasetResults, mlRiskThreshold);

modelDir = fullfile(pwd, 'ml_model');
if ~exist(modelDir,'dir'); mkdir(modelDir); end

save(fullfile(modelDir, 'ml_logistic_model.mat'), 'mlModel', 'trainSummary');

fid = fopen(fullfile(modelDir, 'ml_training_summary.txt'), 'w');
fprintf(fid, 'ML logistic model summary\n');
fprintf(fid, '=========================\n');
fprintf(fid, 'Threshold: %.3f\n', mlRiskThreshold);
fprintf(fid, 'Train rows: %d\n', trainSummary.numTrain);
fprintf(fid, 'Test rows: %d\n', trainSummary.numTest);
fprintf(fid, 'Accuracy: %.4f\n', trainSummary.accuracy);
fprintf(fid, 'Confusion Matrix [TN FP; FN TP]:\n');
fprintf(fid, '[%d %d; %d %d]\n', trainSummary.confMat(1,1), trainSummary.confMat(1,2), trainSummary.confMat(2,1), trainSummary.confMat(2,2));
fprintf(fid, 'Coefficients [Intercept, SNRdB_z, Is64QAM_z, RecentNACKRate_z]:\n');
fprintf(fid, '%.8f ', mlModel.beta);
fprintf(fid, '\n');
fclose(fid);

disp(trainSummary);

%% ------------------------------------------------------------------------
% Phase 3: ML-controlled evaluation sweep
% -------------------------------------------------------------------------
mlCfgList = struct( ...
    'runIdx', {}, ...
    'SNRdB', {}, ...
    'NHARQProcesses', {}, ...
    'rvSeq', {}, ...
    'nTxAnts', {}, ...
    'nRxAnts', {}, ...
    'DelayProfile', {}, ...
    'NumLayers', {}, ...
    'logging', {}, ...
    'dispPlots', {}, ...
    'seed', {}, ...
    'totalNoSlots', {}, ...
    'controlMode', {}, ...
    'fixedModulation', {}, ...
    'datasetMode', {}, ...
    'mlRiskThreshold', {}, ...
    'mlModel', {} );

for k = 1:size(baseIdxList,1)
    cfg = struct();
    cfg.runIdx          = 500000 + k;
    cfg.SNRdB           = SNRdB_vals(baseIdxList(k,1));
    cfg.NHARQProcesses  = NHARQProcesses_vals(baseIdxList(k,2));
    cfg.rvSeq           = rvSeq_vals{baseIdxList(k,3)};
    cfg.nTxAnts         = nTxAnts_vals(baseIdxList(k,4));
    cfg.nRxAnts         = nRxAnts_vals(baseIdxList(k,5));
    cfg.DelayProfile    = DelayProfile_vals{baseIdxList(k,6)};
    cfg.NumLayers       = NumLayers;
    cfg.logging         = 1;
    cfg.dispPlots       = 0;
    cfg.seed            = 300000 + k;
    cfg.totalNoSlots    = totalNoSlots;
    cfg.controlMode     = 'ml';
    cfg.fixedModulation = '16QAM';
    cfg.datasetMode     = false;
    cfg.mlRiskThreshold = mlRiskThreshold;
    cfg.mlModel         = mlModel;
    mlCfgList(end+1) = cfg; %#ok<SAGROW>
end

fprintf('\nPhase 3 ML runs: %d\n', numel(mlCfgList));
mlResults = executeSweep(mlCfgList, numWorkers, 'Phase 3 / ML');

%% ------------------------------------------------------------------------
% Save overall results
% -------------------------------------------------------------------------
save('sweep_results_parallel_ML.mat', 'datasetResults', 'mlResults', 'mlModel', 'trainSummary');

fprintf('\n========================================================\n');
fprintf('  ML workflow complete.\n');
fprintf('  Dataset runs : %d\n', numel(datasetResults));
fprintf('  ML runs      : %d\n', numel(mlResults));
fprintf('  Total time   : %s\n', formatDuration(toc(sweepTimer)));
fprintf('========================================================\n');

%% ========================================================================
% Local functions
% =========================================================================
function results = executeSweep(cfgList, numWorkers, labelText)

totalRuns = numel(cfgList);
results(totalRuns,1) = defaultResultStruct();

fprintf('\n%s\n', repmat('=',1,60));
fprintf('%s\n', labelText);
fprintf('Total runs: %d\n', totalRuns);
fprintf('%s\n', repmat('=',1,60));

dq = parallel.pool.DataQueue;
progressCounter = 0;
phaseStart = tic;
afterEach(dq, @updateProgress);

parfor runIdx = 1:totalRuns
    cfg = cfgList(runIdx);

    try
        if (cfg.nTxAnts < cfg.NumLayers) || (cfg.nRxAnts < cfg.NumLayers)
            r = defaultResultStruct();
            r.runIdx = cfg.runIdx;
            r.ok = false;
            r.skipped = true;
            r.errorMessage = "Invalid antenna/layer combo";
        else
            r = normalizeResultStruct(NR_sim_parallel_ML(cfg));
            r.skipped = false;
            r.errorMessage = "";
        end
    catch ME
        r = defaultResultStruct();
        r.runIdx = cfg.runIdx;
        r.ok = false;
        r.skipped = false;
        r.errorMessage = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
    end

    results(runIdx) = r;
    send(dq, 1);
end

finished = sum(arrayfun(@(x) isfield(x,'ok') && x.ok, results));
skipped  = sum(arrayfun(@(x) isfield(x,'skipped') && x.skipped, results));
failed   = totalRuns - finished - skipped;

fprintf('\n%s complete\n', labelText);
fprintf('Finished: %d | Skipped: %d | Failed: %d | Time: %s\n', ...
    finished, skipped, failed, formatDuration(toc(phaseStart)));

    function updateProgress(~)
        progressCounter = progressCounter + 1;
        if mod(progressCounter,25) == 0 || progressCounter == totalRuns
            fprintf('%s progress: %d / %d completed\n', labelText, progressCounter, totalRuns);
        end
    end
end

function [model, summary] = trainMLModelFromResults(results, threshold)

allTables = {};
for k = 1:numel(results)
    if isfield(results(k),'ok') && results(k).ok && isfield(results(k),'datasetFile') ...
            && strlength(results(k).datasetFile) > 0 && isfile(results(k).datasetFile)
        s = load(results(k).datasetFile, 'datasetTable');
        if isfield(s,'datasetTable') && ~isempty(s.datasetTable)
            allTables{end+1} = s.datasetTable; %#ok<SAGROW>
        end
    end
end

if isempty(allTables)
    failedRuns = results(~[results.ok] & ~[results.skipped]);
    if ~isempty(failedRuns)
        error(['No dataset tables were found because dataset-generation runs failed. ' ...
               'Check the first runtime error above and fix NR_sim_parallel_ML before training.']);
    else
        error('No dataset tables were found. Training cannot continue.');
    end
end

T = vertcat(allTables{:});
T = T(:, {'SNRdB','Is64QAM','RecentNACKRate','NeedRetx'});
T = rmmissing(T);

if height(T) < 10
    error('Dataset is too small for training. Only %d usable rows were found.', height(T));
end

classes = unique(T.NeedRetx);
if numel(classes) < 2
    error('Dataset contains only one class for NeedRetx. Logistic regression requires both 0 and 1 labels.');
end

cv = cvpartition(T.NeedRetx, 'HoldOut', 0.20);
trainTbl = T(training(cv), :);
testTbl  = T(test(cv), :);

Xtrain = trainTbl{:, {'SNRdB','Is64QAM','RecentNACKRate'}};
ytrain = trainTbl.NeedRetx;
Xtest  = testTbl{:, {'SNRdB','Is64QAM','RecentNACKRate'}};
ytest  = testTbl.NeedRetx;

mu = mean(Xtrain,1);
sigma = std(Xtrain,0,1);
sigma(sigma == 0) = 1;

XtrainZ = (Xtrain - mu) ./ sigma;
XtestZ  = (Xtest  - mu) ./ sigma;

trainZ = array2table(XtrainZ, 'VariableNames', {'SNRdB','Is64QAM','RecentNACKRate'});
trainZ.NeedRetx = ytrain;

glm = fitglm(trainZ, 'NeedRetx ~ SNRdB + Is64QAM + RecentNACKRate', 'Distribution','binomial');

beta = glm.Coefficients.Estimate(:);

ptest = sigmoid([ones(size(XtestZ,1),1), XtestZ] * beta);
yhat = double(ptest >= threshold);

cm = confusionmat(ytest, yhat, 'Order', [0 1]);
acc = mean(yhat == ytest);

model = struct();
model.beta = beta;
model.mu = mu;
model.sigma = sigma;
model.featureNames = {'SNRdB','Is64QAM','RecentNACKRate'};
model.threshold = threshold;

summary = struct();
summary.numTrain = numel(ytrain);
summary.numTest = numel(ytest);
summary.accuracy = acc;
summary.confMat = cm;
summary.testMeanProb = mean(ptest);
summary.totalRows = height(T);
end


function r = defaultResultStruct()
r = struct( ...
    'runIdx', 0, ...
    'ok', false, ...
    'skipped', false, ...
    'errorMessage', "", ...
    'logFile', "", ...
    'datasetFile', "", ...
    'controlMode', "", ...
    'SNRdB', nan, ...
    'NHARQProcesses', nan, ...
    'DelayProfile', "", ...
    'ThroughputEfficiencyPct', nan, ...
    'HARQEfficiencyPct', nan, ...
    'RetransmissionRatePct', nan, ...
    'AttemptBLERPct', nan, ...
    'FinalBLERPct', nan, ...
    'AverageThroughputMbps', nan, ...
    'EstimatedLatencyMs', nan, ...
    'Modulation', "" );
end

function r = normalizeResultStruct(in)
r = defaultResultStruct();
fn = fieldnames(r);
for k = 1:numel(fn)
    if isfield(in, fn{k})
        r.(fn{k}) = in.(fn{k});
    end
end
end

function y = sigmoid(x)
y = 1 ./ (1 + exp(-x));
end

function str = formatDuration(seconds)
h = floor(seconds/3600);
m = floor(mod(seconds,3600)/60);
s = mod(seconds,60);
str = sprintf('%02dh:%02dm:%05.2fs', h, m, s);
end
