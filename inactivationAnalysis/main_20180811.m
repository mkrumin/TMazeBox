% this script is for trying things out with inactivations experiments

clear;

list = struct();
iList = 0;

iList = iList + 1;
list(iList).animalName = 'MK027';
list(iList).startDate = '2017-11-08';
list(iList).endDate = '2020-12-01';
list(iList).excludeDate = {'2000-01-01'};
list(iList).excludeSession = {};
% 
iList = iList + 1;
list(iList).animalName = 'JC001';
list(iList).startDate = '2018-04-01';
list(iList).endDate = '2020-12-01';
list(iList).excludeDate = {'2018-06-19'};
list(iList).excludeSession = {}; %date_session_name(expRef)

iList = iList + 1;
list(iList).animalName = 'JC003';
list(iList).startDate = '2018-04-11';
list(iList).endDate = '2020-12-01';
list(iList).excludeDate = {'2018-04-12', '2018-04-19'};
list(iList).excludeSession = {'2018-05-01_2154_JC003', '2018-05-22_1353_JC003'};

% list(iList).startDate = '2018-05-17'; % this is the first day V1 was inactivated
% list(iList).endDate = '2020-12-01';
% list(iList).excludeDate = {'2018-04-12', '2018-04-19'};
% added one session to exclude where V1 was not inactivated yet
list(iList).excludeSession = {'2018-05-01_2154_JC003', '2018-05-22_1353_JC003', '2018-05-17_1353_JC003'};

iList = iList + 1;
list(iList).animalName = 'JC004';
list(iList).startDate = '2018-05-24';
list(iList).endDate = '2020-12-01';
list(iList).excludeDate = {};
list(iList).excludeSession = {'2018-05-30_1117_JC004', '2018-06-05_1019_JC004', '2018-06-06_1043_JC004' }; %date_session_name(expRef)

% iList = iList + 1;
% list(iList).animalName = 'MK031';
% list(iList).startDate = '2018-08-01';
% list(iList).endDate = '2020-12-01';
% list(iList).excludeDate = {};
% list(iList).excludeSession = {}; %date_session_name(expRef)


%%
fprintf('Getting the list of experiments..');
tic
filenames = getFileList(list);
fprintf('.done (%4.2f sec)\n', toc);

%%
fprintf('Loading the data..');
tic
tmazeData = loadTMazeData(filenames);
fprintf('.done (%4.2f sec)\n', toc);

%%
fprintf('Extracting the relevant data..');
tic
res = parseTmazeData(tmazeData);
fprintf('.done (%4.2f sec)\n', toc);

%% Combine all the sessions together
fprintf('Pooling the data..')
tic

contrast = cell2mat({res.contrast}');
outcome = cell2mat({res.outcome}');
behavior = cell2mat({res.behavior}');
finished = cell2mat({res.finished}');
random = cell2mat({res.random}');
optiStim = cell2mat({res.optiStim}');
isV1 = cell2mat({res.isV1}');
isPPC = cell2mat({res.isPPC}');
isUndecided = cell2mat({res.isUndecided}');
thPreStim = cell2mat({res.thPreStim}');
allZ = res(1).z;
allTh = res(1).theta;
nSessions = length(res);
for iSession = 2:nSessions
    allZ = cat(1, allZ, res(iSession).z);
    allTh = cat(1, allTh, res(iSession).theta);
end

fprintf('.done (%4.2f sec)\n', toc);

%% selection of subsets of trials
idxF = finished;
idxRand = random;
idxUnd = isUndecided;

idxNone = ~optiStim(:, 1) & ~optiStim(:,2);
idxLeft = optiStim(:, 1) & ~optiStim(:,2);
idxRight = ~optiStim(:, 1) & optiStim(:,2);
idxBoth = optiStim(:, 1) & optiStim(:,2);

%%
options.figName = reshape([list.animalName], 5, []);
options.figName = reshape([options.figName; repmat(' ', 1, length(list))], 1, []);

options.fitPsycho = true;
nSims = 1000;
nBootSims = 1000;

%% Plotting distribution of thetas at laser stim onset
% all the session and trials are taken into account
figure('name', options.figName, 'Color', [1 1 1]);
histogram(thPreStim);
xlabel('Median-corrected \theta');
ylabel('# trials');
title('\theta at inactivation onset');
xlim([-30 30]);
ax = gca;
ax.FontSize = 16;
box off;
prc = prctile(thPreStim, [25 75]);
hold on;
plot([prc(1) prc(1)], ylim, 'r--', 'LineWidth', 2);
plot([prc(2) prc(2)], ylim, 'r--', 'LineWidth', 2);

%%

clear idx;
idx{1} = idxF & idxRand & idxLeft & isPPC;
idx{2} = idxF & idxRand & idxRight & isPPC;
options.color = {'r', 'b'};
options.groupNames = {'left', 'right'};
options.title = 'All PPC Trials';
options.nSims = nSims;
options.nBootSims = nBootSims;
options.bootParams = [1 1 1 1];
options.testParams = [2 2 1 1];
% parames are: [threshold, slope, guessRate, lapseRate]
% 1 - this parameter will be tested
% 0 - this parameter will be completely unconstrained
% 2 - this parameter will be always constrained between two models 
models = analyzeAndPlot(contrast, behavior, idx, options);
fprintf('Bootstrapping..');
tic;
modBoot = modelBoot(contrast, behavior, idx, models, options);
fprintf('.done (%3.1f sec)\n', toc)
fprintf('Likelihood ratio comparison with Monte-Carlo..');
tic;
modComp = compareModels(contrast, behavior, idx, models, options);
fprintf('.done (%3.1f sec)\n', toc)
printStats(models, modBoot, modComp, options);

%%
clear idx;
idx{1} = idxF & idxRand & idxNone;
idx{2} = idxF & idxRand & idxBoth & isPPC;
options.color = {'k', 'c'};
options.groupNames = {'none', 'both'};
options.title = 'All PPC Trials';
options.nSims = nSims;
options.nBootSims = nBootSims;
options.bootParams = [1 1 1 1];
options.testParams = [2 2 1 1];
% parames are: [threshold, slope, guessRate, lapseRate]
% 1 - this parameter will be tested
% 0 - this parameter will be completely unconstrained
% 2 - this parameter will be always constrained between two models 
models = analyzeAndPlot(contrast, behavior, idx, options);
fprintf('Bootstrapping..');
tic;
modBoot = modelBoot(contrast, behavior, idx, models, options);
fprintf('.done (%3.1f sec)\n', toc)
fprintf('Likelihood ratio comparison with Monte-Carlo..');
tic;
modComp = compareModels(contrast, behavior, idx, models, options);
fprintf('.done (%3.1f sec)\n', toc)
printStats(models, modBoot, modComp, options);


%%
clear idx;
idx{1} = idxF & idxRand & idxLeft & isPPC & isUndecided;
idx{2} = idxF & idxRand & idxRight & isPPC & isUndecided;
options.color = {'r', 'b'};
options.groupNames = {'left', 'right'};
options.title = 'Undecided PPC Trials';
analyzeAndPlot(contrast, behavior, idx, options);

%%
clear idx;
idx{1} = idxF & idxRand & idxLeft & isV1;
idx{2} = idxF & idxRand & idxRight & isV1;
options.color = {'r', 'b'};
options.groupNames = {'left', 'right'};
options.title = 'All V1 Trials';
options.nSims = nSims;
options.nBootSims = 100; nBootSims;
options.bootParams = [1 1 1 1];
options.testParams = [2 2 1 1];
% parames are: [threshold, slope, guessRate, lapseRate]
% 1 - this parameter will be tested
% 0 - this parameter will be completely unconstrained
% 2 - this parameter will be always constrained between two models 
models = analyzeAndPlot(contrast, behavior, idx, options);
fprintf('Bootstrapping..');
tic;
modBoot = modelBoot(contrast, behavior, idx, models, options);
fprintf('.done (%3.1f sec)\n', toc)
fprintf('Likelihood ratio comparison with Monte-Carlo..');
tic;
modComp = compareModels(contrast, behavior, idx, models, options);
fprintf('.done (%3.1f sec)\n', toc)
printStats(models, modBoot, modComp, options);

%%
clear idx;
idx{1} = idxF & idxRand & idxNone;
idx{2} = idxF & idxRand & idxBoth & isV1;
options.color = {'k', 'c'};
options.groupNames = {'none', 'both'};
options.title = 'All V1 Trials';
analyzeAndPlot(contrast, behavior, idx, options);
options.nSims = nSims;
options.nBootSims = nBootSims;
options.bootParams = [1 1 1 1];
options.testParams = [2 1 2 2];
% parames are: [threshold, slope, guessRate, lapseRate]
% 1 - this parameter will be tested
% 0 - this parameter will be completely unconstrained
% 2 - this parameter will be always constrained between two models 
models = analyzeAndPlot(contrast, behavior, idx, options);
fprintf('Bootstrapping..');
tic;
modBoot = modelBoot(contrast, behavior, idx, models, options);
fprintf('.done (%3.1f sec)\n', toc)
fprintf('Likelihood ratio comparison with Monte-Carlo..');
tic;
modComp = compareModels(contrast, behavior, idx, models, options);
fprintf('.done (%3.1f sec)\n', toc)
printStats(models, modBoot, modComp, options);
