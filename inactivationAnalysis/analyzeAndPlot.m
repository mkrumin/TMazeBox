function model = analyzeAndPlot(contrast, behavior, idx, options)

nGroups = length(idx);
model = struct('pars', [], 'LL', [], 'exitflag', [], 'output', []);
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

plotData.cc = cc;
plotData.nr = nr;
plotData.nn = nn;
plotData.pp = pp;
plotData.ci = ci;

plotData.fisher = doFisher(plotData);

%% Fitting psychometric curves

% if options.fitPsycho
%     fprintf('Fitting psychometric curves..');
%     tic
% %     addpath('\\zserver\Code\Psychofit\');
%     nFits = 10;
%     modelType = 'erf_psycho_2gammas';
%     parsMin = [-100, 0, 0, 0];
%     parsMax = [100, 100, 1, 1];
%     xx = [-50:50]';
%     for iGroup = 1:nGroups
%         parsStart = [mean(cc{iGroup}), 10, 0.1, 0.1];
%         [pars L]= mle_fit_psycho([cc{iGroup}, nn{iGroup}, pp{iGroup}]', modelType, ...
%             parsStart, parsMin, parsMax, nFits);
%         yy{iGroup} = erf_psycho_2gammas(pars, xx);
%     end
%     
%    fprintf('.done (%4.2f sec)\n', toc);
% %    rmpath('\\zserver\Code\Psychofit\');
%    plotData.xx = xx;   
%    plotData.yy = yy;
% end
% 

if options.fitPsycho
    fprintf('Fitting psychometric curves..');
    tic
%     addpath('C:\Users\Michael\Documents\MATLAB\Palamedes');
    xx = [-50:50]';
    for iGroup = 1:nGroups
        PF = @PAL_CumulativeNormal;
        searchGrid.alpha = [-50:5:50]; 
        searchGrid.beta = [0:0.01:0.5];
        searchGrid.gamma = [0:0.01:0.4];
        searchGrid.lambda = [0:0.01:0.4];
        paramsFree = [1 1 1 1];
        [pars, LL, exitflag, output] = ...
            PAL_PFML_Fit(cc{iGroup}, nr{iGroup}, nn{iGroup}, searchGrid, paramsFree, PF, ...
            'guessLimits', [0 1], 'lapseLimits', [0 1]);
        model(iGroup).pars = pars;
        model(iGroup).LL = LL;
        model(iGroup).exitflag = exitflag;
        model(iGroup).output = output;
        yy{iGroup} = PAL_CumulativeNormal(pars, xx);
    end
    
   fprintf('.done (%4.2f sec)\n', toc);
%    rmpath('C:\Users\Michael\Documents\MATLAB\Palamedes');
   plotData.xx = xx;   
   plotData.yy = yy;
end


%% Plotting is done here



makeFigure(plotData, options)
