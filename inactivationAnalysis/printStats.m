function printStats(models, modBoot, modComp, options)

hFig = figure('name', [options.figName , options.title], 'Color', [1 1 1]);
hFig.Position = [400 250, 700, 500];
xl = {'Threshold [%]', 'Slope [1/%]', 'Lapse Rate L', 'Lapse Rate R'};
parNames = {'Threshold', 'Slope', 'Lapse L', 'Lapse R'};
nGroups = length(modBoot);
if length(modComp)==1
    % taking care of the case of only a single comparison between two
    % groups
    modComp(1, 2) = modComp(1, 1);
end
for iGroup = 1:nGroups
    idx{iGroup} = modBoot(iGroup).converged;
end
for iPar=1:4
    ax = subplot(2, 2, iPar);
    allVals = [];
    for iGroup = 1:nGroups
        data = modBoot(iGroup).paramSim(idx{iGroup}, iPar);
        hBin(iGroup) = 2*iqr(data)/length(data)^(1/3);
        allVals = cat(1, allVals, data);
    end
    nBins = (max(allVals) - min(allVals))/mean(hBin);
    nBins = max(nBins, 3);
    edges = linspace(min(allVals), max(allVals), nBins + 1);
    for iGroup = 1:nGroups
        histogram(modBoot(iGroup).paramSim(idx{iGroup}, iPar), ...
            edges, 'FaceColor', options.color{iGroup})
        hold on;
    end    %     histogram(modBoot(1).paramSim(:, iPar), nBins, 'FaceColor', options.color{1})
    %     hold on;
    %     histogram(modBoot(2).paramSim(:, iPar), nBins, 'FaceColor', options.color{2})
    if iPar == 1
        leg = legend(options.groupNames);
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
for iGroup = 1:nGroups
    fprintf('\tCondition #%1.0f: %s\n', iGroup, options.groupNames{iGroup});
    fprintf('\t\tThreshold = %5.3f +- %5.3f [%%]\n', models(iGroup).pars(1), modBoot(iGroup).SD(1));
    fprintf('\t\tSlope = %6.4f +- %6.4f [1/%%]\n', models(iGroup).pars(2), modBoot(iGroup).SD(2));
    fprintf('\t\tLapse Left = %6.4f +- %6.4f\n', models(iGroup).pars(3), modBoot(iGroup).SD(3));
    fprintf('\t\tLapse Right = %6.4f +- %6.4f\n', models(iGroup).pars(4), modBoot(iGroup).SD(4));
end

fprintf('\nWilcoxon rank sum test for medians (p-values):\n');
for iGroup = 1:nGroups-1
    for jGroup = iGroup+1:nGroups
    fprintf('\tCondition ''%s'' vs. ''%s'':\n', options.groupNames{iGroup}, options.groupNames{jGroup});
        for iPar = 1:4
            pVal = ranksum(modBoot(iGroup).paramSim(:,iPar), modBoot(jGroup).paramSim(:,iPar));
            fprintf('\t\t%s - %5.3g\n', parNames{iPar}, pVal);
        end
    end
end

fprintf('\nLikelihood ratio model comparison:\n');
fprintf('\tTested parameters: \n')
fprintf('\t\t%s\n', parNames{options.testParams==1})
if sum(options.testParams==0)
    fprintf('\tFree parameters: \n')
    fprintf('\t\t%s\n', parNames{options.testParams==0})
end
if sum(options.testParams==2)
    fprintf('\tConstrained parameters: \n')
    fprintf('\t\t%s\n', parNames{options.testParams==2})
end

nGroups = length(modComp);
for iGroup = 1:nGroups-1
    for jGroup = iGroup+1:nGroups
        fprintf('\tComparing conditions ''%s'' vs. ''%s''\n', ...
            options.groupNames{iGroup}, options.groupNames{jGroup});
        
        fprintf('\tLesser model parameters:\n')
        for iG = [iGroup, jGroup]
            fprintf('\t\tCondition : %s\n', options.groupNames{iG});
            for iPar = 1:4
                fprintf('\t\t\t%s = %6.4f\n', parNames{iPar}, ...
                    modComp(iGroup, jGroup).paramsL(find([iGroup, jGroup] == iG), iPar));
            end
        end
        fprintf('\tFuller model parameters:\n')
        for iG =  [iGroup, jGroup]
            fprintf('\t\tCondition : %s\n', options.groupNames{iG});
            for iPar = 1:4
                fprintf('\t\t\t%s = %6.4f\n', parNames{iPar}, ...
                    modComp(iGroup, jGroup).paramsF(find([iGroup, jGroup] == iG), iPar));
            end
        end
        if modComp(iGroup, jGroup).pTLR > 0
            fprintf('\tP-Value = %5.3d\n', modComp(iGroup, jGroup).pTLR);
        else
            fprintf('\tP-Value < %5.3g\n', 1/length(modComp(iGroup, jGroup).TLRSim));
        end
    end
end
