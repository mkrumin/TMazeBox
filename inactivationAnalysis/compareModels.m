function out = compareModels(contrast, behavior, idx, models, options)

nGroups = length(idx);
alpha = 0.05;

for iGroup = 1:nGroups
    cc{iGroup} = unique(contrast(idx{iGroup}));
    %     cc{iGroup} = cc{iGroup}(~ismember(cc{iGroup}, excludeC));
    nn{iGroup} = nan(size(cc{iGroup}));
    nr{iGroup} = nan(size(cc{iGroup}));
    for iC = 1:length(cc{iGroup})
        nn{iGroup}(iC) = sum(contrast(idx{iGroup}) == cc{iGroup}(iC));
        idxC = idx{iGroup} & contrast == cc{iGroup}(iC);
        nr{iGroup}(iC) = sum(behavior(idxC) == 'R');
    end
    [pp{iGroup}, ci{iGroup}]= binofit(nr{iGroup}, nn{iGroup}, alpha);
end

cLevels = cell2mat(cc)';
nRight = cell2mat(nr)';
nTotal = cell2mat(nn)';
paramValues = cell2mat({models.pars}');
nSims = options.nSims;
optLesser = {'unconstrained', 'constrained', 'constrained'};
optFuller = {'unconstrained', 'unconstrained', 'constrained'};
optVal = options.testParams + 1;
PF = @PAL_CumulativeNormal;
[TLR, pTLR, paramsL, paramsF, TLRSim, converged] = ...
    PAL_PFLR_ModelComparison(cLevels, nRight, nTotal, ...
    paramValues, nSims, PF, 'rangeTries', [20 1 1 1], ...
    'lesserThresholds',  optLesser{optVal(1)}, ...
    'lesserSlopes', optLesser{optVal(2)}, ...
    'lesserGuessrates', optLesser{optVal(3)}, ...
    'lesserLapserates', optLesser{optVal(4)}, ...
    'fullerThresholds', optFuller{optVal(1)}, ...
    'fullerSlopes', optFuller{optVal(2)}, ...
    'fullerGuessrates', optFuller{optVal(3)}, ...
    'fullerLapserates', optFuller{optVal(4)});

out.TLR = TLR;
out.pTLR = pTLR;
out.paramsL = paramsL;
out.paramsF = paramsF;
out.TLRSim = TLRSim;
out.converged = converged;
