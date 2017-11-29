% this script is for trying things out with inactivations experiments

clear;
animalName = 'MK027';
startDate = '2017-11-08';
endDate = '2018-11-29';
excludeDate = {'2000-01-01'};
fitPsycho = true;
alpha = 0.1;
excludeC = NaN;

fprintf('Getting the list of experiments..');
tic
currFolder = cd('\\zserver\Code\Rigging\main');
[ExpRefOld, dateOld] = dat.listExps('MK027');
filenamesOld = dat.expFilePath(ExpRefOld, 'tmaze', 'master');
cd(currFolder);
[ExpRefNew, dateNew] = dat.listExps('MK027');
filenamesNew = dat.expFilePath(ExpRefNew, 'tmaze', 'master');

ExpRef = [ExpRefOld; ExpRefNew];
expDate = [dateOld; dateNew];
filenames = [filenamesOld; filenamesNew];
[expDate, idx] = sort(expDate, 'ascend');
ExpRef = ExpRef(idx);
filenames = filenames(idx);
fprintf('.done (%4.2f sec)\n', toc);

%%
% optogenetic exps didn;t start before that
idxRecent = find(expDate >= datenum(startDate) & ...
    expDate <= datenum(endDate) & ...
    ~ismember(expDate, datenum(excludeDate)));

%%
fprintf('Loading the data..');
tic
tmazeData = struct('EXP', [], 'SESSION', []);
iData = 0;
for iExp = 1:length(idxRecent)
    data = load(filenames{idxRecent(iExp)});
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

idx = finished & random;
% idx = random;
% idx = finished;

idxNone = idx & ~optiStim(:, 1) & ~optiStim(:,2);
idxLeft = idx & optiStim(:, 1) & ~optiStim(:,2);
idxRight = idx & ~optiStim(:, 1) & optiStim(:,2);
idxBoth = idx & optiStim(:, 1) & optiStim(:,2);

idx = {idxNone; idxLeft; idxRight; idxBoth};
if fitPsycho
    LineStyle = {'.k'; '.r'; '.b'; '.m'};
else
    LineStyle = {'.-k'; '.-r'; '.-b'; '.-m'};
end
pcLineStyle = {'k'; 'r'; 'b'; 'm'};
pcLineWidth = 3;
groupName = {'none'; 'left'; 'right'; 'both'};
nGroups = length(idx);
groups2plot = 1:nGroups;
% groups2plot = 2:3;

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


%% Plotting is done here

fprintf('Plotting..')
tic
figure

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
ax.XTick = unique([cc{:}]);
ax.YTick = [0 0.5 1];
ax.XLabel.String = 'Contrast [%]';
ax.YLabel.String = 'Prob (Going Right)';
ax.Title.String = sprintf('%s, starting from %s', animalName, startDate);
axis square;
legend(allLegends);

% write the number of trials near each data point

for iCurve = groups2plot
    for iPoint = 1:length(nn{iCurve})
        tx = text(cc{iCurve}(iPoint)+1, pp{iCurve}(iPoint), sprintf('%1.0f', nn{iCurve}(iPoint)));
        tx.Color = er(iCurve).Color;
        tx.HorizontalAlignment = 'Left';
        tx.VerticalAlignment = 'Middle';
        tx.FontSize = 10;
        tx.FontWeight = 'bold';
    end
end

drawnow;

fprintf('.done (%4.2f sec)\n', toc);