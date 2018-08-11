
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
tic
for iGroup = 1:nGroups
    fprintf('group %1.0f/%1.0f\n', iGroup, nGroups);
    nFiles = length(files{iGroup});
    for iFile = 1:nFiles
        fprintf('file %1.0f/%1.0f\n', iFile, nFiles);
        try
            load(fullfile(folder, files{iGroup}(iFile).name));
        catch
            warning(sprintf('failed to load %s \n', files{iGroup}(iFile).name))
            continue
        end
        
        res{iGroup}(iFile) = getRhoEV(TM);
        
        resFileName = fullfile(folder, 'SuppFigX_rhoev.mat');
        save(resFileName, 'res');
    end
end

resFileName = fullfile(folder, 'SuppFigX_rhoev.mat');
save(resFileName, 'res');

toc

%% load data

folder = 'G:\DATA\';
resFileName = fullfile(folder, 'SuppFigX_rhoev.mat');
load(resFileName);

%% define some parameters
resFull = res;
resFull{end+1} = cell2mat(res);

%% plot rasters
axesLims = [-1 1];
histEdges = axesLims(1):0.05:axesLims(2);
binC = (histEdges(1:end-1) + histEdges(2:end))/2;
baseLineVal = axesLims(1);
histScaling = 2;
fontSize = 14;
groups2plot = 4;

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
    
    figName = sprintf('%s, all data', groupName{iGroup});
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
    plot(axesLims, axesLims, 'k:');
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
    plot(axesLims, axesLims, 'k:');
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
    plot(axesLims, axesLims, 'k:');
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
    plot(axesLims, axesLims, 'k:');
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
    plot(axesLims, axesLims, 'k:');
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
    plot(axesLims, axesLims, 'k:');
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

axesLims = [-1 1];
histEdges = axesLims(1):0.05:axesLims(2);
binC = (histEdges(1:end-1) + histEdges(2:end))/2;
baseLineVal = axesLims(1);
histScaling = 2;
fontSize = 14;
groups2plot = 4;
nGroups = length(res);

figName = sprintf('%s, all data', groupName{iGroup});
hFig = figure('Name', figName, 'Position', [300 500 1200 800]);
nRows = 2;
nColumns = 3;

for iGroup = 1:nGroups
    
    for iSession = 1:length(res{iGroup})
        rhoZThD = res{iGroup}(iSession).rhoZThD;
        rhoZTh = res{iGroup}(iSession).rhoZTh;
        rhoZD = res{iGroup}(iSession).rhoZD;
        evZThD = res{iGroup}(iSession).evZThD_TbyT;
        evZTh = res{iGroup}(iSession).evZTh_TbyT;
        evZD = res{iGroup}(iSession).evZD_TbyT;
        %     evZThD = res{iGroup}.evZThD_SbyS;
        %     evZTh = res{iGroup}.evZTh_SbyS;
        %     evZD = res{iGroup}.evZD_SbyS;
        
        idx = ~isnan(rhoZTh);
        nSelectedCells = sum(idx);
        
        ax = subplot(nRows, nColumns, 1);
%         xMean = mean(rhoZTh(idx));
%         xNeg = std(rhoZTh(idx));
%         xPos = xNeg;
%         yMean = mean(rhoZD(idx));
%         yNeg = std(rhoZD(idx));
%         yPos = yNeg;
        xPr = prctile(rhoZTh(idx), [25 50 75]);
        xData = xPr(2);
        xNeg = xPr(2) - xPr(1);
        xPos = xPr(3) - xPr(2);
        yPr = prctile(rhoZD(idx), [25 50 75]);
        yData = yPr(2);
        yNeg = yPr(2) - yPr(1);
        yPos = yPr(3) - yPr(2);
        er = errorbar(xData, yData, yNeg, yPos, xNeg, xPos, 'o');
        er.MarkerSize = 6;
        er.MarkerFaceColor = er.MarkerEdgeColor;
        er.LineStyle = ':';
        er.LineWidth = 0.25;
        hold on;
        plot(axesLims, axesLims, 'k:');
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
        p = plot(mean(rhoZThD(idx)), mean(rhoZTh(idx)), 'o');
        hold on;
        plot(axesLims, axesLims, 'k:');
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
        p = plot(mean(rhoZThD(idx)), mean(rhoZD(idx)), 'o');
        hold on;
        plot(axesLims, axesLims, 'k:');
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
        plot(axesLims, axesLims, 'k:');
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
        plot(axesLims, axesLims, 'k:');
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
        plot(axesLims, axesLims, 'k:');
        axis equal
        xlim(axesLims)
        ylim(axesLims)
        plot([0 0], axesLims, 'k:', axesLims, [0 0], 'k:');
        xlabel('EV_{z,d,\theta}');
        ylabel('EV_{z,d}')
        ax.FontSize = fontSize;
        box off;
    end
end

%% plot densities with histograms

dr = 0.005;
smoothing = 3;
edges{1} = -1:dr:1;
edges{2} = edges{1};
axesLims = [-1 1];
drHist = 0.02;
edgesHist = -1:drHist:1;
doLogScale = false;

% centres =

cMin = Inf;
cMax = -Inf;
for iGroup = groups2plot
    rhoZThD = cell2mat({resFull{iGroup}.rhoZThD}');
    rhoZTh = cell2mat({resFull{iGroup}.rhoZTh}');
    rhoZD = cell2mat({resFull{iGroup}.rhoZD}');
    zPeak = cell2mat({resFull{iGroup}.zPeak}');
    
    nSessions = length(resFull{iGroup});
    rocArea= [];
    thetaStd = [];
    for iSession = 1:nSessions
        nCells = numel(resFull{iGroup}(iSession).zPeak);
        rocTmp = nan(nCells, 1);
        stdTmp = nan(nCells, 1);
        for iCell = 1:nCells
            iPlane = resFull{iGroup}(iSession).planeIdx(iCell);
            rocTmp(iCell, 1) = interp1(resFull{iGroup}(iSession).zROC{iPlane}, ...
                resFull{iGroup}(iSession).thROC{iPlane}, resFull{iGroup}(iSession).zPeak(iCell));
            stdTmp(iCell, 1) = interp1(resFull{iGroup}(iSession).zROC{iPlane}, ...
                resFull{iGroup}(iSession).thStd{iPlane}, resFull{iGroup}(iSession).zPeak(iCell));
        end
        rocArea = cat(1, rocArea, rocTmp);
        thetaStd = cat(1, thetaStd, stdTmp);
    end
    
    for iZBin = 1:nZBins
        
        %         idx = zPeak>=zBinEdges(iZBin) & zPeak < zBinEdges(iZBin+1);
        %         idx = idx & rocArea<rocThr & thetaStd > stdThr;
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
end

for iGroup = groups2plot
    rhoZThD = cell2mat({resFull{iGroup}.rhoZThD}');
    rhoZTh = cell2mat({resFull{iGroup}.rhoZTh}');
    rhoZD = cell2mat({resFull{iGroup}.rhoZD}');
    zPeak = cell2mat({resFull{iGroup}.zPeak}');
    
    nSessions = length(resFull{iGroup});
    %     rocArea= [];
    %     thetaStd = [];
    %     for iSession = 1:nSessions
    %         nCells = numel(resFull{iGroup}(iSession).zPeak);
    %         rocTmp = nan(nCells, 1);
    %         stdTmp = nan(nCells, 1);
    %         for iCell = 1:nCells
    %             iPlane = resFull{iGroup}(iSession).planeIdx(iCell);
    %             rocTmp(iCell, 1) = interp1(resFull{iGroup}(iSession).zROC{iPlane}, ...
    %                 resFull{iGroup}(iSession).thROC{iPlane}, resFull{iGroup}(iSession).zPeak(iCell));
    %             stdTmp(iCell, 1) = interp1(resFull{iGroup}(iSession).zROC{iPlane}, ...
    %                 resFull{iGroup}(iSession).thStd{iPlane}, resFull{iGroup}(iSession).zPeak(iCell));
    %         end
    %         rocArea = cat(1, rocArea, rocTmp);
    %         thetaStd = cat(1, thetaStd, stdTmp);
    %     end
    %
    for iZBin = 1:nZBins
        
        figName = sprintf('%s, z = [%1.0f %1.0f], AuROC < %4.2f, std_\\theta > %4.2f', groupName{iGroup}, zBinEdges(iZBin), zBinEdges(iZBin+1), rocThr, stdThr);
        figure('Name', figName, 'Position', [300 300 1200 750]);
        cm = colormap('bone');
        colormap(flipud(cm));
        
        %         idx = zPeak>=zBinEdges(iZBin) & zPeak < zBinEdges(iZBin+1);
        %         idx = idx & rocArea<rocThr & thetaStd > stdThr;
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
        
    end
end

