function res = modelBoot(contrast, behavior, idx, models, options)

nGroups = length(idx);
res = struct('SD', [], 'paramSim', [], 'LLSim', [], 'converged', []);
% alpha = 0.05;

for iGroup = 1:nGroups
    cc{iGroup} = unique(contrast(idx{iGroup}));
    %     cc{iGroup} = cc{iGroup}(~ismember(cc{iGroup}, excludeC));
    nn{iGroup} = nan(size(cc{iGroup}));
    %     nr{iGroup} = nan(size(cc{iGroup}));
    for iC = 1:length(cc{iGroup})
        nn{iGroup}(iC) = sum(contrast(idx{iGroup}) == cc{iGroup}(iC));
        idxC = idx{iGroup} & contrast == cc{iGroup}(iC);
        nr{iGroup}(iC) = sum(behavior(idxC) == 'R');
    end
    %     [pp{iGroup}, ci{iGroup}]= binofit(nr{iGroup}, nn{iGroup}, alpha);
end


%% Bootstrapping psychometric curves

for iGroup = 1:nGroups
    PF = @PAL_CumulativeNormal;
    paramsVal = models(iGroup).pars;
    paramsFree = options.bootParams;
    [SD, paramSim, LLSim, converged] = ...
        PAL_PFML_BootstrapParametric(cc{iGroup}(:)', nn{iGroup}(:)', ...
        paramsVal, paramsFree, options.nBootSims, PF, ...
        'guessLimits', [0 1], 'lapseLimits', [0 1]);
    res(iGroup).SD = SD;
    res(iGroup).paramSim = paramSim;
    res(iGroup).LLSim = LLSim;
    res(iGroup).converged = converged;
end


