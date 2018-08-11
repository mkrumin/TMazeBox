function makeFigure(data, options)

fprintf('Plotting..')
tic
hFig = figure('Name', options.figName', 'Position', [650   350   600   500]);
hFig.Color = [1 1 1];

pcLineWidth = 3;
nGroups = length(data.cc);
if options.fitPsycho
    for iGroup = 1:nGroups
        plot(data.xx, data.yy{iGroup}, options.color{iGroup}, 'LineWidth', pcLineWidth);
        hold on;
    end
end

allLegends = cell(0);
% er = nan(length(groups2plot), 1);
for iGroup = 1:nGroups
    er(iGroup) = errorbar(data.cc{iGroup}, data.pp{iGroup}, ...
        data.ci{iGroup}(:,1)-data.pp{iGroup}, ...
        data.ci{iGroup}(:,2)-data.pp{iGroup}, ['.', options.color{iGroup}]);
    hold on;
    er(iGroup).MarkerSize = 30;
    ax = gca;
%     iLegend = find(groups2plot == iGroup);
    allLegends{iGroup} = sprintf('%s (%2.0f trials)', options.groupNames{iGroup}, sum(data.nn{iGroup}));
end

nContrasts = length(data.cc{1});
for iContrast = 1:nContrasts
    if data.fisher.p(iContrast)<0.001
        text(data.cc{1}(iContrast), 0.8, '\Delta', 'FontSize', 18, 'HorizontalAlignment', 'Center');
    elseif data.fisher.p(iContrast)<0.01
        text(data.cc{1}(iContrast), 0.8, '*', 'FontSize', 18, 'HorizontalAlignment', 'Center');
    elseif data.fisher.p(iContrast)<0.05
        text(data.cc{1}(iContrast), 0.8, '+', 'FontSize', 18, 'HorizontalAlignment', 'Center');
    else
        % do nothing
    end
    
end
box off

xlim([-50 50]);
ylim([0 1]);
plot([0 0], ylim, 'k:')
plot(xlim, [0.5 0.5], 'k:')
ax.Color = hFig.Color;
ax.XTick = unique(cell2mat(data.cc(:)));
ax.YTick = [0 0.5 1];
ax.XLabel.String = 'Contrast [%]';
ax.YLabel.String = 'Prob (Going Right)';
ax.FontSize = 16;
ax.Title.String = options.title;
axis square;
legend(allLegends, 'FontSize', 16, 'Location', 'southeast', 'box', 'off');

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

