% this script is for trying things out with inactivations experiments

clear;

list = struct();
iList = 0;

% iList = iList + 1;
% list(iList).animalName = 'MK027';
% list(iList).startDate = '2017-11-08';
% list(iList).endDate = '2020-12-01';
% list(iList).excludeDate = {'2000-01-01'};
% list(iList).excludeSession = {};

% iList = iList + 1;
% list(iList).animalName = 'JC003';
% list(iList).startDate = '2018-04-01';
% list(iList).endDate = '2020-12-01';
% list(iList).excludeDate = {'2018-04-12', '2018-04-19'};
% list(iList).excludeSession = {'2018-05-01_2154_JC003'};

% iList = iList + 1;
% list(iList).animalName = 'JC001';
% list(iList).startDate = '2018-04-01';
% list(iList).endDate = '2020-12-01';
% list(iList).excludeDate = {};
% list(iList).excludeSession = {};

iList = iList + 1;
list(iList).animalName = 'JC004';
list(iList).startDate = '2018-05-01';
list(iList).endDate = '2020-12-01';
list(iList).excludeDate = {};
list(iList).excludeSession = {};

fitPsycho = true;
fitSmartPC = false;
alpha = 0.05;
excludeC = NaN;
% groups2plot = 1:3;
groups2plot = 2:3;
% groups2plot = 2:4;
% groups2plot = [1, 4];
% groups2plot = 1:4;

%%
fprintf('Getting the list of experiments..');
filenames = {};
for iList = 1:length(list)
    tic
    currFolder = cd('\\zserver\Code\Rigging\main');
    [ExpRefOld, dateOld] = dat.listExps(list(iList).animalName);
    filenamesOld = dat.expFilePath(ExpRefOld, 'tmaze', 'master');
    cd(currFolder);
    [ExpRefNew, dateNew] = dat.listExps(list(iList).animalName);
    filenamesNew = dat.expFilePath(ExpRefNew, 'tmaze', 'master');
    
    ExpRef = [ExpRefOld; ExpRefNew];
    expDate = [dateOld; dateNew];
    fn = [filenamesOld; filenamesNew];
    [expDate, idx] = sort(expDate, 'ascend');
    ExpRef = ExpRef(idx);
    fn = fn(idx);
    fprintf('.done (%4.2f sec)\n', toc);
    
    % optogenetic exps didn't start before that
    idxRecent = expDate >= datenum(list(iList).startDate) & expDate <= datenum(list(iList).endDate);
    % exclude certain dates dates
    idxRecent = idxRecent & ~ismember(expDate, datenum(list(iList).excludeDate));
    % exclude certain sessions
    idxRecent = idxRecent & ~ismember(ExpRef, list(iList).excludeSession);
    
    idxRecent = find(idxRecent);
    
    filenames = cat(1, filenames, fn(idxRecent));
    
end
%%
fprintf('Loading the data..');
tic
tmazeData = struct('EXP', [], 'SESSION', []);
iData = 0;
for iExp = 1:length(filenames)
    data = load(filenames{iExp});
    % only take experiments with inactivations
    if data.EXP.optiStim
        iData = iData + 1;
        tmazeData(iData) = data;
    end
end
fprintf('.done (%4.2f sec)\n', toc);

%% Combine all the sessions together
fprintf('Pooling the data..')
tic
res = struct;
nSessions = length(tmazeData);
for iSession = 1:nSessions
    SESSION = tmazeData(iSession).SESSION;
    EXP = tmazeData(iSession).EXP;
    nTrials = length(SESSION.allTrials);
    res(iSession).contrast = nan(nTrials, 1);
    res(iSession).outcome = '';
    res(iSession).behavior = '';
    res(iSession).finished = nan(nTrials, 1);
    res(iSession).random = nan(nTrials, 1);
    % HACK assuming two locations
    res(iSession).optiStim = nan(nTrials, 2);
    res(iSession).z = cell(nTrials, 1);
    res(iSession).theta  = cell(nTrials, 1);
    
    for iTrial = 1:nTrials
        res(iSession).contrast(iTrial) = SESSION.allTrials(iTrial).info.contrast;
        side = -1+2*isequal(SESSION.allTrials(iTrial).info.stimulus, 'RIGHT');
        res(iSession).contrast(iTrial) = res(iSession).contrast(iTrial) * side;
        res(iSession).outcome(iTrial) = SESSION.allTrials(iTrial).info.outcome(1);
        res(iSession).behavior(iTrial) = res(iSession).outcome(iTrial);
        if res(iSession).outcome(iTrial) == 'C'
            res(iSession).behavior(iTrial) = SESSION.allTrials(iTrial).info.stimulus(1);
        elseif res(iSession).outcome(iTrial) == 'W'
            res(iSession).behavior(iTrial) = char('R'+'L'-SESSION.allTrials(iTrial).info.stimulus(1));
        end
        res(iSession).finished(iTrial) = ismember(res(iSession).behavior(iTrial), {'R', 'L'});
        if isequal(EXP.stimType, 'BAITED')
            res(iSession).random(iTrial) = iTrial == 1 || res(iSession).outcome(iTrial-1)=='C';
        else
            res(iSession).random(iTrial) = true; %assuming 'RANDOM'
        end
        res(iSession).optiStim(iTrial, :) = SESSION.allTrials(iTrial).info.optiStim.laserPower;
        zInd = find(ismember(SESSION.allTrials(iTrial).pospars, 'Z'));
        thInd = find(ismember(SESSION.allTrials(iTrial).pospars, 'theta'));
        res(iSession).z{iTrial} = -SESSION.allTrials(iTrial).posdata(:, zInd);
        res(iSession).theta{iTrial} = SESSION.allTrials(iTrial).posdata(:, thInd)*180/pi;
    end
    res(iSession).behavior = res(iSession).behavior(:);
    res(iSession).outcome = res(iSession).outcome(:);
end
res = res(:);

% combine all the sessions

contrast = cell2mat({res.contrast}');
outcome = cell2mat({res.outcome}');
behavior = cell2mat({res.behavior}');
finished = cell2mat({res.finished}');
random = cell2mat({res.random}');
optiStim = cell2mat({res.optiStim}');
allZ = res(1).z;
allTh = res(1).theta;
for iSession = 2:nSessions
    allZ = cat(1, allZ, res(iSession).z);
    allTh = cat(1, allTh, res(iSession).theta);
end

idx = finished & random;
% idx = random;
% idx = finished;

idxNone = idx & ~optiStim(:, 1) & ~optiStim(:,2);
idxLeft = idx & optiStim(:, 1) & ~optiStim(:,2);
idxRight = idx & ~optiStim(:, 1) & optiStim(:,2);
idxBoth = idx & optiStim(:, 1) & optiStim(:,2);

idx = {idxNone; idxLeft; idxRight; idxBoth};
if fitPsycho
    LineStyle = {'.k'; '.r'; '.b'; '.c'};
else
    LineStyle = {'.-k'; '.-r'; '.-b'; '.-c'};
end
pcLineStyle = {'k'; 'r'; 'b'; 'c'};
pcLineWidth = 2;
groupName = {'none'; 'left'; 'right'; 'both'};

for iGroup = groups2plot
    cc{iGroup} = unique(contrast(idx{iGroup}));
    cc{iGroup} = cc{iGroup}(~ismember(cc{iGroup}, excludeC));
    nn{iGroup} = nan(size(cc{iGroup}));
    nr{iGroup} = nan(size(cc{iGroup}));
    for iC = 1:length(cc{iGroup})
        nn{iGroup}(iC) = sum(contrast(idx{iGroup}) == cc{iGroup}(iC));
        idxC = idx{iGroup} & contrast == cc{iGroup}(iC);
        nr{iGroup}(iC) = sum(behavior(idxC) == 'R');
    end
    [pp{iGroup}, ci{iGroup}]= binofit(nr{iGroup}, nn{iGroup}, alpha);
end

fprintf('.done (%4.2f sec)\n', toc);

%% Fitting psychometric curves

if fitPsycho
    fprintf('Fitting psychometric curves..');
    tic
    addpath('\\zserver\Code\Psychofit\');
    nFits = 10;
    modelType = 'erf_psycho_2gammas';
    parsMin = [-100, 0, 0, 0];
    parsMax = [100, 100, 1, 1];
    xx = [-50:50]';
    for iGroup = groups2plot
        parsStart = [mean(cc{iGroup}), 10, 0.1, 0.1];
        [pars L]= mle_fit_psycho([cc{iGroup}, nn{iGroup}, pp{iGroup}]', modelType, ...
            parsStart, parsMin, parsMax, nFits);
        yy{iGroup} = erf_psycho_2gammas(pars, xx);
    end
    
    fprintf('.done (%4.2f sec)\n', toc);
    rmpath('\\zserver\Code\Psychofit\');
end

%% Fitting psychometric curves with corrected biases
if fitSmartPC
    nTrials = cellfun(@length, {res.contrast}', 'UniformOutput', 0);
    iSession = num2cell([1:length(nTrials)]');
    iSession = cellfun(@(x, y) ones(x, 1)*y, nTrials, iSession, 'UniformOutput', 0);
    iSession = cell2mat(iSession);
    
    allData = struct;
    allData.contrast = contrast;
    allData.outcome = outcome;
    allData.behavior = behavior;
    allData.finished = finished;
    allData.random = random;
    allData.optiStim = optiStim;
    allData.iSession = iSession;
    
    dataPC = fitDebiasedPC(allData);
end


%% Plotting is done here

fprintf('Plotting..')
tic
hFig = figure;
hFig.Color = [1 1 1];

if fitPsycho
    for iGroup = groups2plot
        plot(xx, yy{iGroup}, pcLineStyle{iGroup}, 'LineWidth', pcLineWidth);
        hold on;
    end
end

allLegends = cell(0);
% er = nan(length(groups2plot), 1);
for iGroup = groups2plot
    er(iGroup) = errorbar(cc{iGroup}, pp{iGroup}, ...
        ci{iGroup}(:,1)-pp{iGroup}, ci{iGroup}(:,2)-pp{iGroup}, LineStyle{iGroup});
    er(iGroup).MarkerSize = 30;
    ax = gca;
    iLegend = find(groups2plot == iGroup);
    allLegends{iLegend} = sprintf('%s (%2.0f)', groupName{iGroup}, sum(nn{iGroup}));
end
box off

xlim([-50 50]);
ylim([0 1]);
plot([0 0], ylim, 'k:')
plot(xlim, [0.5 0.5], 'k:')
ax.Color = hFig.Color;
ax.XTick = unique([cc{:}]);
ax.YTick = [0 0.5 1];
ax.XLabel.String = 'Contrast [%]';
ax.YLabel.String = 'Prob (Going Right)';
% ax.Title.String = sprintf('%s, starting from %s', animalName, startDate);
axis square;
legend(allLegends);

% write the number of trials near each data point

% for iCurve = groups2plot
%     for iPoint = 1:length(nn{iCurve})
%         tx = text(cc{iCurve}(iPoint)+1, pp{iCurve}(iPoint), sprintf('%1.0f', nn{iCurve}(iPoint)));
%         tx.Color = er(iCurve).Color;
%         tx.HorizontalAlignment = 'Left';
%         tx.VerticalAlignment = 'Middle';
%         tx.FontSize = 10;
%         tx.FontWeight = 'bold';
%     end
% end
% 
% ccAll = unique([cc{:}]);
% nnAll = zeros(size(ccAll));
% for iC = 1:length(ccAll)
%     for iGroup = groups2plot
%         nnAll(iC) = nnAll(iC) + sum(nn{iGroup}(cc{iGroup} == ccAll(iC)));
%     end
%     tx = text(ccAll(iC), 0, sprintf('%1.0f', nnAll(iC)));
%     tx.Color = [0 0 0];
%     tx.HorizontalAlignment = 'Center';
%     tx.VerticalAlignment = 'Bottom';
%     tx.FontSize = 12;
%     tx.FontWeight = 'bold';
%     
% end

drawnow;

fprintf('.done (%4.2f sec)\n', toc);

%% plotting the trajectories

prc = [25 50 75];
zLims = [5 60];
cc = unique(abs(contrast));
idxFR = finished & random;
idxCorrect = outcome == 'C';
idxR = behavior == 'R';
idxL = behavior == 'L';
idxStimL = optiStim(:, 1) > 0 & optiStim(:, 2) == 0;
idxStimR = optiStim(:, 2) > 0 & optiStim(:, 1) == 0;
idxStimBoth = optiStim(:, 1) > 0 & optiStim(:, 2) > 0;
idxStimNone = optiStim(:, 1) == 0 & optiStim(:, 2) == 0;
[zAxis, thMatrix] = binZTh(allZ, allTh, 30);

figure('Name', list(1).animalName)

nRows = length(cc)-1;
for iC = 2:length(cc)
    idxC = (abs(contrast) == cc(iC));
    
    subplot(nRows, 4, 1+(iC-2)*4)
    idxBlue = find(idxC & idxCorrect & idxR & idxStimL & idxFR);
    idxRed = find(idxC & idxCorrect & idxL & idxStimL & idxFR);
    
    thLeft = prctile(thMatrix(:, idxRed), prc, 2);
    thRight = prctile(thMatrix(:, idxBlue), prc, 2);
    plot(thLeft(:, 2), zAxis, 'r', 'LineWidth', 2)
    hold on;
    plot(thLeft(:, [1, 3]), zAxis, 'r:', 'LineWidth', 1);
    plot(thRight(:, 2), zAxis, 'b', 'LineWidth', 2);
    plot(thRight(:, [1, 3]), zAxis, 'b:', 'LineWidth', 1);
    axis equal tight
    ylim(zLims);
    title(sprintf('Stim L (%1.0f; %1.0f)', length(idxRed), length(idxBlue)));
    ylabel(sprintf('\\pm%1.0f %%', cc(iC)));
    
    subplot(nRows, 4, 2+(iC-2)*4)
    idxBlue = find(idxC & idxCorrect & idxR & idxStimR & idxFR);
    idxRed = find(idxC & idxCorrect & idxL & idxStimR & idxFR);
    
    thLeft = prctile(thMatrix(:, idxRed), prc, 2);
    thRight = prctile(thMatrix(:, idxBlue), prc, 2);
    plot(thLeft(:, 2), zAxis, 'r', 'LineWidth', 2)
    hold on;
    plot(thLeft(:, [1, 3]), zAxis, 'r:', 'LineWidth', 1);
    plot(thRight(:, 2), zAxis, 'b', 'LineWidth', 2);
    plot(thRight(:, [1, 3]), zAxis, 'b:', 'LineWidth', 1);
    axis equal tight
    ylim(zLims);
    title(sprintf('Stim R (%1.0f; %1.0f)', length(idxRed), length(idxBlue)));
    
    subplot(nRows, 4, 3+(iC-2)*4)
    idxBlue = find(idxC & idxCorrect & idxR & idxStimBoth & idxFR);
    idxRed = find(idxC & idxCorrect & idxL & idxStimBoth & idxFR);
    
    thLeft = prctile(thMatrix(:, idxRed), prc, 2);
    thRight = prctile(thMatrix(:, idxBlue), prc, 2);
    plot(thLeft(:, 2), zAxis, 'r', 'LineWidth', 2)
    hold on;
    plot(thLeft(:, [1, 3]), zAxis, 'r:', 'LineWidth', 1);
    plot(thRight(:, 2), zAxis, 'b', 'LineWidth', 2);
    plot(thRight(:, [1, 3]), zAxis, 'b:', 'LineWidth', 1);
    axis equal tight
    ylim(zLims);
    title(sprintf('Stim Both (%1.0f; %1.0f)', length(idxRed), length(idxBlue)));
    
    subplot(nRows, 4, 4+(iC-2)*4)
    idxBlue = find(idxC & idxCorrect & idxR & idxStimNone & idxFR);
    idxRed = find(idxC & idxCorrect & idxL & idxStimNone & idxFR);
    
    thLeft = prctile(thMatrix(:, idxRed), prc, 2);
    thRight = prctile(thMatrix(:, idxBlue), prc, 2);
    plot(thLeft(:, 2), zAxis, 'r', 'LineWidth', 2)
    hold on;
    plot(thLeft(:, [1, 3]), zAxis, 'r:', 'LineWidth', 1);
    plot(thRight(:, 2), zAxis, 'b', 'LineWidth', 2);
    plot(thRight(:, [1, 3]), zAxis, 'b:', 'LineWidth', 1);
    axis equal tight
    ylim(zLims);
    title(sprintf('Stim None (%1.0f; %1.0f)', length(idxRed), length(idxBlue)));
    
end

%% plotting the trajectories - all in one

prc = [25 50 75];
zLims = [5 60];
cc = unique(abs(contrast));
idxFR = finished & random;
idxCorrect = outcome == 'C';
idxR = behavior == 'R';
idxL = behavior == 'L';
idxStimL = optiStim(:, 1) > 0 & optiStim(:, 2) == 0;
idxStimR = optiStim(:, 2) > 0 & optiStim(:, 1) == 0;
idxStimBoth = optiStim(:, 1) > 0 & optiStim(:, 2) > 0;
idxStimNone = optiStim(:, 1) == 0 & optiStim(:, 2) == 0;
[zAxis, thMatrix] = binZTh(allZ, allTh, 30);

figure('Name', list(1).animalName)

nRows = 1;
nC = length(cc);
shade = [nC-1:-1:0]'/nC;
for iC = 1:length(cc)
    idxC = (abs(contrast) == cc(iC));
    
    subplot(nRows, 4, 1)
    if cc(iC) == 0
        idxBlue = find(idxC & idxR & idxStimL & idxFR);
        idxRed = find(idxC & idxL & idxStimL & idxFR);
    else
        idxBlue = find(idxC & idxCorrect & idxR & idxStimL & idxFR);
        idxRed = find(idxC & idxCorrect & idxL & idxStimL & idxFR);
    end
    
    thLeft = prctile(thMatrix(:, idxRed), prc, 2);
    thRight = prctile(thMatrix(:, idxBlue), prc, 2);
    plot(thLeft(:, 2), zAxis, 'r', 'LineWidth', 2)
    hold on;
    plot(thLeft(:, [1, 3]), zAxis, 'r:', 'LineWidth', 1);
    plot(thRight(:, 2), zAxis, 'b', 'LineWidth', 2);
    plot(thRight(:, [1, 3]), zAxis, 'b:', 'LineWidth', 1);
    axis equal tight
    ylim(zLims);
    title(sprintf('Stim L (%1.0f; %1.0f)', length(idxRed), length(idxBlue)));
    ylabel(sprintf('\\pm%1.0f %%', cc(iC)));
    
    subplot(nRows, 4, 2+(iC-2)*4)
    idxBlue = find(idxC & idxCorrect & idxR & idxStimR & idxFR);
    idxRed = find(idxC & idxCorrect & idxL & idxStimR & idxFR);
    
    thLeft = prctile(thMatrix(:, idxRed), prc, 2);
    thRight = prctile(thMatrix(:, idxBlue), prc, 2);
    plot(thLeft(:, 2), zAxis, 'r', 'LineWidth', 2)
    hold on;
    plot(thLeft(:, [1, 3]), zAxis, 'r:', 'LineWidth', 1);
    plot(thRight(:, 2), zAxis, 'b', 'LineWidth', 2);
    plot(thRight(:, [1, 3]), zAxis, 'b:', 'LineWidth', 1);
    axis equal tight
    ylim(zLims);
    title(sprintf('Stim R (%1.0f; %1.0f)', length(idxRed), length(idxBlue)));
    
    subplot(nRows, 4, 3+(iC-2)*4)
    idxBlue = find(idxC & idxCorrect & idxR & idxStimBoth & idxFR);
    idxRed = find(idxC & idxCorrect & idxL & idxStimBoth & idxFR);
    
    thLeft = prctile(thMatrix(:, idxRed), prc, 2);
    thRight = prctile(thMatrix(:, idxBlue), prc, 2);
    plot(thLeft(:, 2), zAxis, 'r', 'LineWidth', 2)
    hold on;
    plot(thLeft(:, [1, 3]), zAxis, 'r:', 'LineWidth', 1);
    plot(thRight(:, 2), zAxis, 'b', 'LineWidth', 2);
    plot(thRight(:, [1, 3]), zAxis, 'b:', 'LineWidth', 1);
    axis equal tight
    ylim(zLims);
    title(sprintf('Stim Both (%1.0f; %1.0f)', length(idxRed), length(idxBlue)));
    
    subplot(nRows, 4, 4+(iC-2)*4)
    idxBlue = find(idxC & idxCorrect & idxR & idxStimNone & idxFR);
    idxRed = find(idxC & idxCorrect & idxL & idxStimNone & idxFR);
    
    thLeft = prctile(thMatrix(:, idxRed), prc, 2);
    thRight = prctile(thMatrix(:, idxBlue), prc, 2);
    plot(thLeft(:, 2), zAxis, 'r', 'LineWidth', 2)
    hold on;
    plot(thLeft(:, [1, 3]), zAxis, 'r:', 'LineWidth', 1);
    plot(thRight(:, 2), zAxis, 'b', 'LineWidth', 2);
    plot(thRight(:, [1, 3]), zAxis, 'b:', 'LineWidth', 1);
    axis equal tight
    ylim(zLims);
    title(sprintf('Stim None (%1.0f; %1.0f)', length(idxRed), length(idxBlue)));
    
end

%%

cc = unique(abs(contrast));
idxF = finished;
idxRan = random;
idxCorrect = outcome == 'C';
idxR = behavior == 'R';
idxL = behavior == 'L';
idxStimL = optiStim(:, 1) > 0 & optiStim(:, 2) == 0;
idxStimR = optiStim(:, 2) > 0 & optiStim(:, 1) == 0;
idxStimBoth = optiStim(:, 1) > 0 & optiStim(:, 2) > 0;
idxStimNone = optiStim(:, 1) == 0 & optiStim(:, 2) == 0;
trialDur = nan(size(allZ));
for iTrial = 1:length(allZ)
    trialDur(iTrial) = sum(allZ{iTrial}>5)/60;
end

%%
figure('Name', list(1).animalName)
nFinBoth = sum(idxRan & idxF & idxStimBoth);
nBoth = sum(idxRan & idxStimBoth);
[pBoth, pciBoth] = binofit(nFinBoth, nBoth, alpha);
nFinNone = sum(idxRan & idxF & idxStimNone); 
nNone = sum(idxRan & idxStimNone);
[pNone, pciNone] = binofit(nFinNone, nNone, alpha);
nFinRight = sum(idxRan & idxF & idxStimR); 
nRight = sum(idxRan & idxStimR);
[pRight, pciRight] = binofit(nFinRight, nRight, alpha);
nFinLeft = sum(idxRan & idxF & idxStimL); 
nLeft = sum(idxRan & idxStimL);
[pLeft, pciLeft] = binofit(nFinLeft, nLeft, alpha);

% histogram(trialDur(idxRan & idxStimNone))
probs = [pNone, pRight, pLeft, pBoth];
pciLow = probs - [pciNone(1), pciRight(1), pciLeft(1), pciBoth(1)];
pciHigh = [pciNone(2), pciRight(2), pciLeft(2), pciBoth(2)] - probs;
hBar = errorbar([1:4], probs, pciLow, pciHigh, 'o');
ax = hBar.Parent;
ax.XTick = 1:4;
ax.XTickLabel = {'none', 'right', 'left', 'both'};
xlim([0.5, 4.5]);
title('Probability of finishing a random trial');
xlabel('Inactivation condition');
ylabel('Prob to finish with 95% c.i.');
box off

figure('Name', list(1).animalName)
meanBoth = mean(trialDur(idxRan & idxF & idxStimBoth));
stdBoth = std(trialDur(idxRan & idxF & idxStimBoth))/sqrt(nFinBoth);
meanNone = mean(trialDur(idxRan & idxF & idxStimNone));
stdNone = std(trialDur(idxRan & idxF & idxStimNone))/sqrt(nFinNone);
meanRight = mean(trialDur(idxRan & idxF & idxStimR));
stdRight = std(trialDur(idxRan & idxF & idxStimR))/sqrt(nFinRight);
meanLeft = mean(trialDur(idxRan & idxF & idxStimL));
stdLeft = std(trialDur(idxRan & idxF & idxStimL))/sqrt(nFinLeft);

durs = [meanNone, meanRight, meanLeft, meanBoth];
errs = [stdNone, stdRight, stdLeft, stdBoth];
hBar = errorbar([1:4], durs, errs, 'o');
ax = hBar.Parent;
ax.XTick = 1:4;
ax.XTickLabel = {'none', 'right', 'left', 'both'};
xlim([0.5, 4.5]);
title('Mean trial duration');
xlabel('Inactivation condition');
ylabel('trial duration [s] \pm s.e.m.');
box off

figure('Name', list(1).animalName)
[prcBoth] = prctile(trialDur(idxRan & idxF & idxStimBoth), [25 50 75]);
[prcNone] = prctile(trialDur(idxRan & idxF & idxStimNone), [25 50 75]);
[prcRight] = prctile(trialDur(idxRan & idxF & idxStimR), [25 50 75]);
[prcLeft] = prctile(trialDur(idxRan & idxF & idxStimL), [25 50 75]);

durs = [prcNone(2), prcRight(2), prcLeft(2), prcBoth(2)];
durLow = durs - [prcNone(1), prcRight(1), prcLeft(1), prcBoth(1)];
durHigh = [prcNone(1), prcRight(1), prcLeft(1), prcBoth(1)] - durs;
hBar = errorbar([1:4], durs, durLow, durHigh, 'o');
ax = hBar.Parent;
ax.XTick = 1:4;
ax.XTickLabel = {'none', 'right', 'left', 'both'};
xlim([0.5, 4.5]);
title('Median \pm quartiles trial duration');
xlabel('Inactivation condition');
ylabel('trial duration [s]');
box off
