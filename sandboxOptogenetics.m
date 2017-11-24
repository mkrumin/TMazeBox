% this script is for trying things out with inactivations experiments

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

%%
% optogenetic exps didn;t start before that
idxRecent = find(expDate >= datenum('2017-11-08'));


%%
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
    
%%
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

%% combine all the sessions

contrast = cell2mat({res.contrast}');
outcome = cell2mat({res.outcome}');
behavior = cell2mat({res.behavior}');
finished = cell2mat({res.finished}');
random = cell2mat({res.random}');
optiStim = cell2mat({res.optiStim}');

idx = finished & random;
% idx = finished;

idxNone = idx & ~optiStim(:, 1) & ~optiStim(:,2);
idxLeft = idx & optiStim(:, 1) & ~optiStim(:,2);
idxRight = idx & ~optiStim(:, 1) & optiStim(:,2);
idxBoth = idx & optiStim(:, 1) & optiStim(:,2);

%%
alpha = 0.05;
excludeC = [-25, 25];

ccNone = unique(contrast(idxNone));
ccNone = ccNone(~ismember(ccNone, excludeC));
nnNone = nan(size(ccNone));
nrNone = nan(size(ccNone));
for iC = 1:length(ccNone)
    nnNone(iC) = sum(contrast(idxNone) == ccNone(iC));
    idxC = idxNone & contrast == ccNone(iC);
    nrNone(iC) = sum(behavior(idxC) == 'R');
end
[ppNone, ciNone]= binofit(nrNone, nnNone, alpha);

ccLeft = unique(contrast(idxLeft));
ccLeft = ccLeft(~ismember(ccLeft, excludeC));
nnLeft = nan(size(ccLeft));
nrLeft = nan(size(ccLeft));
for iC = 1:length(ccLeft)
    nnLeft(iC) = sum(contrast(idxLeft) == ccLeft(iC));
    idxC = idxLeft & contrast == ccLeft(iC);
    nrLeft(iC) = sum(behavior(idxC) == 'R');
end
[ppLeft, ciLeft]= binofit(nrLeft, nnLeft, alpha);

ccRight = unique(contrast(idxRight));
ccRight = ccRight(~ismember(ccRight, excludeC));
nnRight = nan(size(ccRight));
nrRight = nan(size(ccRight));
for iC = 1:length(ccRight)
    nnRight(iC) = sum(contrast(idxRight) == ccRight(iC));
    idxC = idxRight & contrast == ccRight(iC);
    nrRight(iC) = sum(behavior(idxC) == 'R');
end
[ppRight, ciRight]= binofit(nrRight, nnRight, alpha);

ccBoth = unique(contrast(idxBoth));
ccBoth = ccBoth(~ismember(ccBoth, excludeC));
nnBoth = nan(size(ccBoth));
nrBoth = nan(size(ccBoth));
for iC = 1:length(ccBoth)
    nnBoth(iC) = sum(contrast(idxBoth) == ccBoth(iC));
    idxC = idxBoth & contrast == ccBoth(iC);
    nrBoth(iC) = sum(behavior(idxC) == 'R');
end
[ppBoth, ciBoth]= binofit(nrBoth, nnBoth, alpha);

%%

figure
errorbar(ccNone, ppNone, ciNone(:,1)-ppNone, ciNone(:,2)-ppNone, 'o-k')
ax = gca;
hold on;
errorbar(ccLeft, ppLeft, ciLeft(:,1)-ppLeft, ciLeft(:,2)-ppLeft, 'o-r')
errorbar(ccRight, ppRight, ciRight(:,1)-ppRight, ciRight(:,2)-ppRight, 'o-b')
% errorbar(ccBoth, ppBoth, ciBoth(:,1)-ppBoth, ciBoth(:,2)-ppBoth, 'o-')
box off

xlim([-50 50]);
ylim([0 1]);
plot([0 0], ylim, 'k:')
plot(xlim, [0.5 0.5], 'k:')
ax.XTick = unique([ccNone; ccLeft; ccRight; ccBoth]);
ax.YTick = [0 0.5 1];
legend('none', 'left', 'right', 'both');
