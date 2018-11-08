function printStats(models, modBoot, modComp, options)

hFig = figure('name', [options.figName , options.title], 'Color', [1 1 1]);
hFig.Position = [400 250, 700, 500];
xl = {'Threshold [%]', 'Slope [1/%]', 'Lapse Rate L', 'Lapse Rate R'};
parNames = {'Threshold', 'Slope', 'Lapse L', 'Lapse R'};
nBins = 30;
for iPar=1:4
    ax = subplot(2, 2, iPar);
    idx1 = modBoot(1).converged;
    idx2 = modBoot(2).converged;
    allVals = cat(1, modBoot(1).paramSim(idx1, iPar), modBoot(2).paramSim(idx2, iPar));
    edges = linspace(min(allVals), max(allVals), nBins + 1);
    histogram(modBoot(1).paramSim(idx1, iPar), edges, 'FaceColor', options.color{1})
    hold on;
    histogram(modBoot(2).paramSim(idx2, iPar), edges, 'FaceColor', options.color{2})
    %     histogram(modBoot(1).paramSim(:, iPar), nBins, 'FaceColor', options.color{1})
    %     hold on;
    %     histogram(modBoot(2).paramSim(:, iPar), nBins, 'FaceColor', options.color{2})
    if iPar == 1
        leg = legend(options.groupNames{1}, options.groupNames{2});
        leg.Box = 'off';
        leg.Location = 'best';
        title(options.title)
    end
    xlabel(xl{iPar});
    box off;
    ax.FontSize = 16;
    axis tight
end

fprintf('Animals included: %s\n', options.figName);
fprintf('Trials included: %s\n', options.title);
fprintf('Estimated parameters, with s.e.m. from bootstrapping:\n');
for iGroup = 1:2
    fprintf('\tCondition #%1.0f: %s\n', iGroup, options.groupNames{iGroup});
    fprintf('\t\tThreshold = %5.3f +- %5.3f [%%]\n', models(iGroup).pars(1), modBoot(iGroup).SD(1));
    fprintf('\t\tSlope = %6.4f +- %6.4f [1/%%]\n', models(iGroup).pars(2), modBoot(iGroup).SD(2));
    fprintf('\t\tLapse Left = %6.4f +- %6.4f\n', models(iGroup).pars(3), modBoot(iGroup).SD(3));
    fprintf('\t\tLapse Right = %6.4f +- %6.4f\n', models(iGroup).pars(4), modBoot(iGroup).SD(4));
end

fprintf('\nWilcoxon rank sum test for medians (p-values):\n');
for iPar = 1:4
    pVal = ranksum(modBoot(1).paramSim(:,iPar), modBoot(2).paramSim(:,iPar));
    fprintf('\t%s - %d\n', parNames{iPar}, pVal);
end

fprintf('\nLikelihood ratio model comparison:\n');
fprintf('Tested parameters: \n')
fprintf('\t%s\n', parNames{options.testParams==1})
if sum(options.testParams==0)
    fprintf('Free parameters: \n')
    fprintf('\t%s\n', parNames{options.testParams==0})
end
if sum(options.testParams==2)
    fprintf('Constrained parameters: \n')
    fprintf('\t%s\n', parNames{options.testParams==2})
end
fprintf('Lesser model parameters:\n')
for iGroup = 1:2
    fprintf('\tCondition : %s\n', options.groupNames{iGroup});
    for iPar = 1:4
        fprintf('\t\t%s = %6.4f\n', parNames{iPar}, modComp.paramsL(iGroup, iPar));
    end
end
fprintf('Fuller model parameters:\n')
for iGroup = 1:2
    fprintf('\tCondition : %s\n', options.groupNames{iGroup});
    for iPar = 1:4
        fprintf('\t\t%s = %6.4f\n', parNames{iPar}, modComp.paramsF(iGroup, iPar));
    end
end
if modComp.pTLR > 0
    fprintf('P-Value = %4.2d\n', modComp.pTLR);
else
    fprintf('P-Value < %4.2d\n', 1/length(modComp.TLRSim));
end