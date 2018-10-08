
folder = 'G:\DATA\';

aavDatasets = {...
    '2014-08-02_2023_MK012_TMwExtras.mat'; ...
    '2014-08-04_1827_MK012_TMwExtras.mat'; ...
    '2014-08-05_2228_MK012_TMwExtras.mat'; ...
    '2014-08-08_2251_MK012_TMwExtras.mat'; ...
    '2014-08-11_2133_MK012_TMwExtras.mat'; ...
    '2014-08-13_2023_MK012_TMwExtras.mat'; ...
    '2014-08-15_1931_MK012_TMwExtras.mat'; ...
    '2014-08-02_2203_MK014_TMwExtras.mat'; ...
    '2014-08-05_1937_MK014_TMwExtras.mat'; ...
    };

vglutDatasets = {...
    '2017-06-12_1420_JL005_TMwExtras.mat'; ...
    '2017-07-15_1708_JL008_TMwExtras.mat'; ...
    '2017-07-27_1433_JL008_TMwExtras.mat'; ...
    '2017-08-12_1056_JL008_TMwExtras.mat'; ...
    '2017-09-18_1707_JL008_TMwExtras.mat'; ...
    '2017-09-23_1539_JL008_TMwExtras.mat'; ...
    };

tripleDatasets = {...
    '2015-07-03_2127_MK020_TMwExtras.mat'; ...
    '2015-07-07_2127_MK020_TMwExtras.mat'; ...
    '2015-07-30_2010_MK020_TMwExtras.mat'; ...
    '2015-07-31_1851_MK020_TMwExtras.mat'; ...
    '2015-08-02_1709_MK020_TMwExtras.mat'; ...
    '2015-08-05_1930_MK020_TMwExtras.mat'; ...
    '2015-08-07_1821_MK020_TMwExtras.mat'; ...
    '2015-07-03_1632_MK022_TMwExtras.mat'; ...
    '2015-08-02_2051_MK022_TMwExtras.mat'; ...
    '2015-06-19_2155_MK023_TMwExtras.mat'; ...
    '2015-07-18_154_MK023_TMwExtras.mat';
    };

allPPC = [aavDatasets; tripleDatasets; vglutDatasets];

allFiles = dir(sprintf('%s*_TMwExtras.mat', folder));
files{1} = allFiles(ismember({allFiles.name}, aavDatasets));
files{2} = allFiles(ismember({allFiles.name}, tripleDatasets));
files{3} = allFiles(ismember({allFiles.name}, vglutDatasets));

% idx = ismember({allFiles.name}, allPPC);
% % idx = ismember({allFiles.name}, '2017-09-23_1539_JL008_TMwExtras.mat');
% % idx = ismember({allFiles.name}, '2015-08-02_1709_MK020_TMwExtras.mat');
% idx = ismember({allFiles.name}, '2014-08-05_1937_MK014_TMwExtras.mat');
%
% files = mat2cell(allFiles(idx), ones(size(allFiles(idx))), 1);

nGroups = length(files);

groupName = {'AAV', '3tg', 'VGlut', 'All mice'};
% for iGroup = 1:nGroups
%     groupName{iGroup} = files{iGroup}(1).name(1:21);
% end
groupName{nGroups+1} = 'All mice';
groups2plot = 1:length(groupName);

%% analysis

% tic
% for iGroup = 1:nGroups
%     fprintf('group %1.0f/%1.0f\n', iGroup, nGroups);
%     nFiles = length(files{iGroup});
%     for iFile = 1:nFiles
%         fprintf('file %1.0f/%1.0f\n', iFile, nFiles);
%         try
%             load(fullfile(folder, files{iGroup}(iFile).name));
%         catch
%             warning(sprintf('failed to load %s \n', files{iGroup}(iFile).name))
%             continue
%         end
%         
%         res{iGroup}(iFile) = getRhoEV_wAuROC(TM);
%         
%         resFileName = fullfile(folder, 'SuppFigX_rhoev_wAuROC.mat');
%         save(resFileName, 'res');
%     end
% end
% 
% resFileName = fullfile(folder, 'SuppFigX_rhoev_wAuROC.mat');
% save(resFileName, 'res');
% 
% toc

%% load data

folder = 'G:\DATA\';
resFileName = fullfile(folder, 'SuppFigX_rhoev_wAuROC.mat');
whatData = 'ROC < 0.95 samples';
% resFileName = fullfile(folder, 'SuppFigX_rhoev.mat');
% whatData = 'all samples';
load(resFileName);

%% define some parameters
resFull = res;
resFull{end+1} = cell2mat(res);

%% plot scatters
axesLims = [-.5 1];
histEdges = axesLims(1):0.05:axesLims(2);
binC = (histEdges(1:end-1) + histEdges(2:end))/2;
baseLineVal = axesLims(1);
histScaling = 2;
fontSize = 14;
groups2plot = 3;
diagThickness = 2;

for iGroup = groups2plot
    
    rhoZThD = cell2mat({resFull{iGroup}.rhoZThD}');
    rhoZTh = cell2mat({resFull{iGroup}.rhoZTh}');
    rhoZD = cell2mat({resFull{iGroup}.rhoZD}');
    evZThD = cell2mat({resFull{iGroup}.evZThD_TbyT}');
    evZTh = cell2mat({resFull{iGroup}.evZTh_TbyT}');
    evZD = cell2mat({resFull{iGroup}.evZD_TbyT}');
    %     evZThD = cell2mat({resFull{iGroup}.evZThD_SbyS}');
    %     evZTh = cell2mat({resFull{iGroup}.evZTh_SbyS}');
    %     evZD = cell2mat({resFull{iGroup}.evZD_SbyS}');
    
    figName = sprintf('%s, %s', groupName{iGroup}, whatData);
    hFig = figure('Name', figName, 'Position', [300 500 1200 800]);
    nRows = 2;
    nColumns = 3;
    
    %         idx = zPeak>=zBinEdges(iZBin) & zPeak < zBinEdges(iZBin+1);
    %         idx = idx & rocArea<rocThr & thetaStd > stdThr;
    idx = ~isnan(rhoZTh);
    nSelectedCells = sum(idx);
    
    ax = subplot(nRows, nColumns, 1);
    p = plot(rhoZTh(idx), rhoZD(idx), '.');
    hold on;
    plot(axesLims, axesLims, 'k:', 'LineWidth', diagThickness);
    axis equal
    xlim(axesLims)
    ylim(axesLims)
    plot([0 0], axesLims, 'k:', axesLims, [0 0], 'k:');
    xlabel('\rho_{z,\theta}');
    ylabel('\rho_{z,d}')
    box off;
    text(axesLims(1)+0.1, 0.9, sprintf('nCells = %1.0f', nSelectedCells), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 16)
    ax.FontSize = fontSize;
    
    ax = subplot(nRows, nColumns, 2);
    p = plot(rhoZThD(idx), rhoZTh(idx), '.');
    hold on;
    plot(axesLims, axesLims, 'k:', 'LineWidth', diagThickness);
    axis equal
    xlim(axesLims)
    ylim(axesLims)
    plot([0 0], axesLims, 'k:', axesLims, [0 0], 'k:');
    xlabel('\rho_{z,d,\theta}')
    ylabel('\rho_{z,\theta}');
    title(figName);
    ax.FontSize = fontSize;
    box off;
    
    ax = subplot(nRows, nColumns, 3);
    p = plot(rhoZThD(idx), rhoZD(idx), '.');
    hold on;
    plot(axesLims, axesLims, 'k:', 'LineWidth', diagThickness);
    axis equal
    xlim(axesLims)
    ylim(axesLims)
    plot([0 0], axesLims, 'k:', axesLims, [0 0], 'k:');
    xlabel('\rho_{z,d,\theta}');
    ylabel('\rho_{z,d}')
    ax.FontSize = fontSize;
    box off;
    
    idx = ~isnan(evZTh);
    nSelectedCells = sum(idx);
    
    ax = subplot(nRows, nColumns, 4);
    plot(evZTh(idx), evZD(idx), '.');
    hold on;
    plot(axesLims, axesLims, 'k:', 'LineWidth', diagThickness);
    axis equal
    xlim(axesLims)
    ylim(axesLims)
    plot([0 0], axesLims, 'k:', axesLims, [0 0], 'k:');
    xlabel('EV_{z,\theta}');
    ylabel('EV_{z,d}')
    box off;
    text(axesLims(1)+0.1, 0.9, sprintf('nCells = %1.0f', nSelectedCells), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 16)
    ax.FontSize = fontSize;
    
    ax = subplot(nRows, nColumns, 5);
    plot(evZThD(idx), evZTh(idx), '.');
    hold on;
    plot(axesLims, axesLims, 'k:', 'LineWidth', diagThickness);
    axis equal
    xlim(axesLims)
    ylim(axesLims)
    plot([0 0], axesLims, 'k:', axesLims, [0 0], 'k:');
    xlabel('EV_{z,d,\theta}')
    ylabel('EV_{z,\theta}');
    ax.FontSize = fontSize;
    box off;
    
    ax = subplot(nRows, nColumns, 6);
    plot(evZThD(idx), evZD(idx), '.');
    hold on;
    plot(axesLims, axesLims, 'k:', 'LineWidth', diagThickness);
    axis equal
    xlim(axesLims)
    ylim(axesLims)
    plot([0 0], axesLims, 'k:', axesLims, [0 0], 'k:');
    xlabel('EV_{z,d,\theta}');
    ylabel('EV_{z,d}')
    ax.FontSize = fontSize;
    box off;
end

% plot marginal histograms for all the axes
hCh = hFig.Children;
for iCh=1:length(hCh)
    if isequal(hCh(iCh).Type, 'axes')
        h = hCh(iCh).Children;
        axes(hCh(iCh));
        for i=1:length(h)
            if isequal(h(i).Type, 'line') && isequal(h(i).Marker, '.')
                [N] = histcounts(h(i).XData, histEdges);
                %                 N = N/sum(N)/mean(diff(histEdges))/histScaling;
                N = N/max(N)/histScaling;
                b = bar(binC, baseLineVal+N);
                b.FaceColor = 'none';
                b.BarWidth = 1;
                b.ShowBaseLine = 'off';
                b.BaseValue = baseLineVal;
                [N] = histcounts(h(i).YData, histEdges);
                %                 N = N/sum(N)/mean(diff(histEdges))/histScaling;
                N = N/max(N)/histScaling;
                b = barh(binC, baseLineVal+N);
                b.FaceColor = 'none';
                b.BarWidth = 1;
                b.ShowBaseLine = 'off';
                b.BaseValue = baseLineVal;
                hCh(iCh).Clipping = 'off';
                
            end
        end
    end
end

%% plot session-by-session summary

axesLims = [-.05 .5];
ticks = [0:0.1:0.5];
histEdges = axesLims(1):0.05:axesLims(2);
binC = (histEdges(1:end-1) + histEdges(2:end))/2;
baseLineVal = axesLims(1);
histScaling = 2;
fontSize = 14;
groups2plot = 4;
nGroups = length(resFull);
plotType = 'median'; % {'mean', 'median', 'meanstd', 'medianprc'}

figName = sprintf('%s, %s', groupName{iGroup}, whatData);
hFig = figure('Name', figName, 'Position', [300 200 1200 800]);
nRows = 2;
nColumns = 3;

for iGroup = groups2plot
    
    for iSession = 1:length(resFull{iGroup})
        rhoZThD = resFull{iGroup}(iSession).rhoZThD;
        rhoZTh = resFull{iGroup}(iSession).rhoZTh;
        rhoZD = resFull{iGroup}(iSession).rhoZD;
        evZThD = resFull{iGroup}(iSession).evZThD_TbyT;
        evZTh = resFull{iGroup}(iSession).evZTh_TbyT;
        evZD = resFull{iGroup}(iSession).evZD_TbyT;
        %     evZThD = res{iGroup}.evZThD_SbyS;
        %     evZTh = res{iGroup}.evZTh_SbyS;
        %     evZD = res{iGroup}.evZD_SbyS;
        
        idx = ~isnan(rhoZTh);
        nSelectedCells = sum(idx);
        
        ax = subplot(nRows, nColumns, 1);
        plotSbySSummary(ax, rhoZTh(idx), rhoZD(idx), plotType, axesLims)
        xlabel('\rho_{z,\theta}');
        ylabel('\rho_{z,d}')
        box off;
        %         text(axesLims(1)+0.1, 0.9, sprintf('nCells = %1.0f', nSelectedCells), ...
        %             'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 16)
        ax.FontSize = fontSize;
        ax.XTick = ticks;
        
        ax = subplot(nRows, nColumns, 2);
        plotSbySSummary(ax, rhoZThD(idx), rhoZTh(idx), plotType, axesLims)
        xlabel('\rho_{z,d,\theta}')
        ylabel('\rho_{z,\theta}');
        title(figName);
        ax.FontSize = fontSize;
        ax.XTick = ticks;
        box off;
        
        ax = subplot(nRows, nColumns, 3);
        plotSbySSummary(ax, rhoZThD(idx), rhoZD(idx), plotType, axesLims)
        xlabel('\rho_{z,d,\theta}');
        ylabel('\rho_{z,d}')
        ax.FontSize = fontSize;
        ax.XTick = ticks;
        box off;
        
        idx = ~isnan(evZTh);
        nSelectedCells = sum(idx);
        
        ax = subplot(nRows, nColumns, 4);
        plotSbySSummary(ax, evZTh(idx), evZD(idx), plotType, axesLims)
        xlabel('EV_{z,\theta}');
        ylabel('EV_{z,d}')
        box off;
        %         text(axesLims(1)+0.1, 0.9, sprintf('nCells = %1.0f', nSelectedCells), ...
        %             'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 16)
        ax.FontSize = fontSize;
        ax.XTick = ticks;
        
        ax = subplot(nRows, nColumns, 5);
        plotSbySSummary(ax, evZThD(idx), evZTh(idx), plotType, axesLims)
        xlabel('EV_{z,d,\theta}')
        ylabel('EV_{z,\theta}');
        ax.FontSize = fontSize;
        ax.XTick = ticks;
        box off;
        
        ax = subplot(nRows, nColumns, 6);
        plotSbySSummary(ax, evZThD(idx), evZD(idx), plotType, axesLims)
        xlabel('EV_{z,d,\theta}');
        ylabel('EV_{z,d}')
        ax.FontSize = fontSize;
        ax.XTick = ticks;
        box off;
    end
end

%% plot densities with histograms

dr = 0.005;
smoothing = 3;
edges{1} = -1:dr:1;
edges{2} = edges{1};
axesLims = [-0 1];
drHist = 0.02;
edgesHist = -1:drHist:1;
doLogScale = false;

cMin = Inf;
cMax = -Inf;
for iGroup = groups2plot
    rhoZThD = cell2mat({resFull{iGroup}.rhoZThD}');
    rhoZTh = cell2mat({resFull{iGroup}.rhoZTh}');
    rhoZD = cell2mat({resFull{iGroup}.rhoZD}');
    
    nSessions = length(resFull{iGroup});
    
    idx = ~isnan(rhoZTh);
    
    [N, C] = hist3([rhoZD(idx), rhoZTh(idx)], 'Edges', edges);
    N = imgaussfilt(N, smoothing);
    N = N/sum(N(:))/dr^2;
    if doLogScale
        N = log10(N+1);
    end
    cMin = min(cMin, min(N(:)));
    cMax = max(cMax, max(N(:)));
    [N, C] = hist3([rhoZThD(idx), rhoZTh(idx)], 'Edges', edges);
    N = imgaussfilt(N, smoothing);
    N = N/sum(N(:))/dr^2;
    if doLogScale
        N = log10(N+1);
    end
    cMin = min(cMin, min(N(:)));
    cMax = max(cMax, max(N(:)));
    [N, C] = hist3([rhoZD(idx), rhoZThD(idx)], 'Edges', edges);
    N = imgaussfilt(N, smoothing);
    N = N/sum(N(:))/dr^2;
    if doLogScale
        N = log10(N+1);
    end
    cMin = min(cMin, min(N(:)));
    cMax = max(cMax, max(N(:)));
    
end

for iGroup = groups2plot
    rhoZThD = cell2mat({resFull{iGroup}.rhoZThD}');
    rhoZTh = cell2mat({resFull{iGroup}.rhoZTh}');
    rhoZD = cell2mat({resFull{iGroup}.rhoZD}');
    
    nSessions = length(resFull{iGroup});
    
    figName = sprintf('%s', groupName{iGroup});
    figure('Name', figName, 'Position', [300 300 1200 750]);
    cm = colormap('bone');
    colormap(flipud(cm));
    
    idx = ~isnan(rhoZTh);
    nSelectedCells = sum(idx);
    
    ax = subplot(2, 3, 1);
    [N, C] = hist3([rhoZD(idx), rhoZTh(idx)], 'Edges', edges);
    N = imgaussfilt(N, smoothing);
    N = N/sum(N(:))/dr^2;
    if doLogScale
        N = log10(N+1);
    end
    imagesc(C{2}(1:end-1), C{1}(1:end-1), N(1:end-1, 1:end-1));
    hold on;
    plot(axesLims, axesLims, 'r--', 'LineWidth', 1);
    axis equal xy
    xlim(axesLims)
    ylim(axesLims)
    plot([0 0], axesLims, 'k:', axesLims, [0 0], 'k:');
    xlabel('\rho_{z,\theta}');
    ylabel('\rho_{z,d}')
    %         caxis([cMin cMax]);
    box off;
    ax.XTick  = [0, 0.5, 1];
    ax.YTick  = [0, 0.5, 1];
    ax.FontSize = 14;
    text(0.1, 0.9, sprintf('nCells = %1.0f', nSelectedCells), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 16)
    
    ax = subplot(2, 3, 4);
    [counts, centers] = hist(rhoZTh(idx)-rhoZD(idx), -1:drHist:1);
    c50 = prctile(rhoZTh(idx)-rhoZD(idx), 50);
    cm = nanmean(rhoZTh(idx)-rhoZD(idx));
    counts = counts/sum(counts(:))/drHist/2;
    %     stairs(centers, log10(counts+1));
    stairs(centers, counts);
    hold on;
    %         plot([c50 c50], ylim, 'r--')
    %         plot([cm cm], ylim, 'g--')
    axis tight
    plot([0 0], ylim, 'k:')
    xlim([-0.5 0.5]);
    %         legend('', 'p50', 'mean')
    xlabel('\rho_{z,\theta} - \rho_{z,d}');
    box off;
    ax.YTick  = [];
    ax.FontSize = 14;
    
    ax = subplot(2, 3, 2);
    [N, C] = hist3([rhoZTh(idx), rhoZThD(idx)], 'Edges', edges);
    N = imgaussfilt(N, smoothing);
    N = N/sum(N(:))/dr^2;
    if doLogScale
        N = log10(N+1);
    end
    imagesc(C{2}(1:end-1), C{1}(1:end-1), N(1:end-1, 1:end-1));
    hold on;
    plot(axesLims, axesLims, 'r--', 'LineWidth', 1);
    axis equal xy
    xlim(axesLims)
    ylim(axesLims)
    plot([0 0], axesLims, 'k:', axesLims, [0 0], 'k:');
    xlabel('\rho_{z,d,\theta}')
    ylabel('\rho_{z,\theta}');
    %         caxis([cMin cMax]);
    box off;
    title(figName);
    ax.XTick  = [0, 0.5, 1];
    ax.YTick  = [0, 0.5, 1];
    ax.FontSize = 14;
    
    ax = subplot(2, 3, 5);
    [counts, centers] = hist(rhoZThD(idx)-rhoZTh(idx), -1:drHist:1);
    c50 = prctile(rhoZTh(idx)-rhoZD(idx), 50);
    cm = nanmean(rhoZTh(idx)-rhoZD(idx));
    counts = counts/sum(counts(:))/drHist/2;
    %     stairs(centers, log10(counts+1));
    stairs(centers, counts);
    hold on;
    %         plot([c50 c50], ylim, 'r--')
    %         plot([cm cm], ylim, 'g--')
    box off;
    ax.YTick  = [];
    ax.FontSize = 14;
    axis tight
    plot([0 0], ylim, 'k:')
    xlim([-0.5 0.5]);
    xlabel('\rho_{z,d,\theta} - \rho_{z,\theta}')
    
    ax = subplot(2, 3, 3);
    [N, C] = hist3([rhoZD(idx), rhoZThD(idx)], 'Edges', edges);
    N = imgaussfilt(N, smoothing);
    N = N/sum(N(:))/dr^2;
    if doLogScale
        N = log10(N+1);
    end
    imagesc(C{2}(1:end-1), C{1}(1:end-1), N(1:end-1, 1:end-1));
    hold on;
    plot(axesLims, axesLims, 'r--', 'LineWidth', 1);
    axis equal xy
    xlim(axesLims)
    ylim(axesLims)
    plot([0 0], axesLims, 'k:', axesLims, [0 0], 'k:');
    xlabel('\rho_{z,d,\theta}');
    ylabel('\rho_{z,d}')
    %         caxis([cMin cMax]);
    box off;
    ax.XTick  = [0, 0.5, 1];
    ax.YTick  = [0, 0.5, 1];
    ax.FontSize = 14;
    
    ax = subplot(2, 3, 6);
    [counts, centers] = hist(rhoZThD(idx)-rhoZD(idx), -1:drHist:1);
    c50 = prctile(rhoZTh(idx)-rhoZD(idx), 50);
    cm = nanmean(rhoZTh(idx)-rhoZD(idx));
    counts = counts/sum(counts(:))/drHist/2;
    %     stairs(centers, log10(counts+1));
    stairs(centers, counts);
    hold on;
    %         plot([c50 c50], ylim, 'r--')
    %         plot([cm cm], ylim, 'g--')
    axis tight
    plot([0 0], ylim, 'k:')
    xlim([-0.5 0.5]);
    xlabel('\rho_{z,d,\theta} - \rho_{z,d}');
    box off;
    ax.YTick  = [];
    ax.FontSize = 14;
    
    %     end
end

%% calculate stats of the results

dr = 0.005;
smoothing = 3;
edges{1} = -1:dr:1;
edges{2} = edges{1};
axesLims = [-0 1];
drHist = 0.02;
edgesHist = -1:drHist:1;

groups2plot = 4;

for iGroup = groups2plot
    rhoZThD = cell2mat({resFull{iGroup}.rhoZThD}');
    rhoZTh = cell2mat({resFull{iGroup}.rhoZTh}');
    rhoZD = cell2mat({resFull{iGroup}.rhoZD}');
    
    idx = ~isnan(rhoZTh);
    
    dRho = rhoZThD(idx)-rhoZTh(idx);
    figure;
    ax = subplot(1, 1, 1);
    histogram(dRho, edgesHist)
    xlim([-0.4 0.4])
    dRhoMedian = median(dRho);
    dRhoMAD = mad(dRho, 1);
    pVal = signtest(dRho);
    %     p = signtest(rhoZThD(idx), rhoZTh(idx))
    xlabel('\rho_{z,d,\theta} - \rho_{z,\theta}')
    hold on;
    plot([0 0], ylim, 'k:')
    plot([dRhoMedian, dRhoMedian], ylim, 'r--')
    title({sprintf('\\Delta\\rho = %4.2d \\pm %4.2d (median\\pmm.a.d.)', dRhoMedian, dRhoMAD);...
        sprintf('p-value (two-sided sign-test) = %4.2d', pVal)});
    ax.FontSize = 14;
    box off
    ax.YTick = [];
    ax.YTickLabel = '';
    
    dRho = rhoZTh(idx)-rhoZD(idx);
    figure;
    ax = subplot(1, 1, 1);
    histogram(dRho, edgesHist)
    xlim([-0.4 0.4])
    dRhoMedian = median(dRho);
    dRhoMAD = mad(dRho, 1);
    pVal = signtest(dRho);
    %     p = signtest(rhoZThD(idx), rhoZTh(idx))
    xlabel('\rho_{z,\theta} - \rho_{z,d}')
    hold on;
    plot([0 0], ylim, 'k:')
    plot([dRhoMedian, dRhoMedian], ylim, 'r--')
    title({sprintf('\\Delta\\rho = %4.2d \\pm %4.2d (median\\pmm.a.d.)', dRhoMedian, dRhoMAD);...
        sprintf('p-value (two-sided sign-test) = %4.2d', pVal)});
    ax.FontSize = 14;
    box off
    ax.YTick = [];
    ax.YTickLabel = '';

    dRho = rhoZThD(idx)-rhoZD(idx);
    figure;
    ax = subplot(1, 1, 1);
    histogram(dRho, edgesHist)
    xlim([-0.4 0.4])
    dRhoMedian = median(dRho);
    dRhoMAD = mad(dRho, 1);
    pVal = signtest(dRho);
    %     p = signtest(rhoZThD(idx), rhoZTh(idx))
    xlabel('\rho_{z,d,\theta} - \rho_{z,d}')
    hold on;
    plot([0 0], ylim, 'k:')
    plot([dRhoMedian, dRhoMedian], ylim, 'r--')
    title({sprintf('\\Delta\\rho = %4.2d \\pm %4.2d (median\\pmm.a.d.)', dRhoMedian, dRhoMAD);...
        sprintf('p-value (two-sided sign-test) = %4.2d', pVal)});
    ax.FontSize = 14;
    box off
    ax.YTick = [];
    ax.YTickLabel = '';

end

%% plot marginal histograms for different genotypes
axesLims = [-0.5 1];
histEdges = axesLims(1):0.05:axesLims(2);
binC = (histEdges(1:end-1) + histEdges(2:end))/2;
baseLineVal = axesLims(1);
histScaling = 2;
fontSize = 12;
groups2plot = 1:4;
lineThickness = 2;

figName = sprintf('By genotype, %s', whatData);
hFig = figure('Name', figName, 'Position', [20 200 1800 400]);
for iGroup = groups2plot
    
    rhoZTh = cell2mat({resFull{iGroup}.rhoZTh}');
    
    nRows = 1;
    nColumns = length(groups2plot);
    
    idx = ~isnan(rhoZTh);
    nSelectedCells = sum(idx);
    
    ax = subplot(nRows, nColumns, find(groups2plot==iGroup));
    hHist = histogram(rhoZTh(idx), histEdges);
    hold on;
    axis square
    xlim(axesLims)
    yLimits = ylim;
    ylim([0 yLimits(2)])
    
    %     ylim(axesLims)
    rhoMedian = median(rhoZTh(idx));
    rhoMAD = mad(rhoZTh(idx), 1);
    plot([0 0], ylim, 'k--', 'LineWidth', 2);
    plot([1 1]*rhoMedian, ylim, 'r--', 'LineWidth', 2);
    xlabel('\rho_{z,\theta}');
    box off;
    text(axesLims(1)+0.1, yLimits(2), sprintf('nCells = %1.0f', nSelectedCells), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 16)
    ax.FontSize = fontSize;
    title({groupName{iGroup}; sprintf('%4.2f\\pm%4.2f (med\\pmmad)', rhoMedian, rhoMAD)})
    ax.YTick = [];
    
end

%% plotting histograms-by-genotype on the same axes
hFig = figure('Name', figName);
clrs = [1 0 0; 0 1 0; 0 0 1];
for iGroup = 1:3
    
    rhoZTh = cell2mat({resFull{iGroup}.rhoZTh}');
    
    idx = ~isnan(rhoZTh);
    nSelectedCells = sum(idx);
    
    %     hHist = histogram(rhoZTh(idx), histEdges, 'Normalization', 'pdf', 'FaceColor', 'none', 'EdgeColor', clrs(iGroup, :));
    hHist = histogram(rhoZTh(idx), histEdges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2);
    hold on;
end
ax = gca;
axis square
xlim(axesLims)
yLimits = ylim;
ylim([0 yLimits(2)])

%     ylim(axesLims)
plot([0 0], ylim, 'k--', 'LineWidth', 2);
xlabel('\rho_{z,\theta}');
box off;
ax.FontSize = fontSize;
ax.YTick = [];
legend(groupName(1:3));

%% comparing medians
for iGroup = 1:3
    r{iGroup} = cell2mat({resFull{iGroup}.rhoZTh}');
    idx = ~isnan(r{iGroup});
    r{iGroup} = r{iGroup}(idx);
end

p = nan(3);
for i=1:3
    for j=i+1:3
        p(i, j) = ranksum(r{i}, r{j});
    end
end

%% plot marginal histograms for different mice

allExpRefs = {resFull{end}.ExpRef}';
animalNames = {'MK012', 'MK014', 'MK020', 'MK022', 'MK023', 'JL005', 'JL008'};
genotype = {'AAV','AAV', '3g', '3g', '3g', 'vGlut','vGlut'};

animalIdx = cellfun(@contains, ...
    repmat(allExpRefs, 1, length(animalNames)),...
    repmat(animalNames, length(allExpRefs), 1));
% add column of ones as last 'animal' - for all animals
animalIdx = cat(2, animalIdx, ones(size(animalIdx, 1), 1));
animalNames = [animalNames, {'All'}];
genotype = [genotype, {'Mix'}];
nAnimals = length(animalNames);

axesLims = [-0.5 1];
histEdges = axesLims(1):0.05:axesLims(2);
binC = (histEdges(1:end-1) + histEdges(2:end))/2;
baseLineVal = axesLims(1);
histScaling = 2;
fontSize = 12;
groups2plot = 1:4;
lineThickness = 2;

figName = sprintf('By mouse, %s', whatData);
hFig = figure('Name', figName, 'Position', [20 200 1800 800]);
for iAnimal = 1:nAnimals
    
    iSession = find(animalIdx(:, iAnimal));
    tmp = resFull{end};
    rhoZTh = cell2mat({tmp(iSession).rhoZTh}');
    
    nRows = 2;
    nColumns = ceil(nAnimals/nRows);
    
    idx = ~isnan(rhoZTh);
    nCells = sum(idx);
    
    ax = subplot(nRows, nColumns, iAnimal);
    hHist = histogram(rhoZTh(idx), histEdges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
    hold on;
    axis square
    xlim(axesLims)
    yLimits = ylim;
    ylim([0 yLimits(2)])
    
    %     ylim(axesLims)
    rhoMedian(iAnimal) = median(rhoZTh(idx));
    rhoMAD(iAnimal) = mad(rhoZTh(idx), 1);
    group{iAnimal} = genotype{iAnimal};
    plot([0 0], ylim, 'k--', 'LineWidth', 2);
    plot([1 1]*rhoMedian(iAnimal), ylim, 'r--', 'LineWidth', 2);
    xlabel('\rho_{z,\theta}');
    box off;
    text(axesLims(2), yLimits(2), sprintf('%1.0f neurons\n%4.2f+-%4.2f', ...
        nCells, rhoMedian(iAnimal), rhoMAD(iAnimal)), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 16)
    ax.FontSize = fontSize;
%     title({animalNames{iAnimal}; sprintf('%4.2f\\pm%4.2f (med\\pmmad)', rhoMedian(iAnimal), rhoMAD(iAnimal))})
    title(animalNames{iAnimal});
    ax.YTick = [];
    
end

anova1(rhoMedian(1:end-1), genotype(1:end-1))