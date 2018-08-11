
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

%%
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
        
        res{iGroup}(iFile) = getFitQuality(TM);
       
        resFileName = fullfile(folder, 'Fig3Groups_EVNoOptZ_1.mat');
        save(resFileName, 'res');
    end
end

resFileName = fullfile(folder, 'Fig3Groups_EVNoOptZ_1.mat');
save(resFileName, 'res');

toc

%% load data

folder = 'G:\DATA\';
% resFileName = fullfile(folder, 'Fig3Groups_RhoNoOptZ_095.mat');
% resFileName = fullfile(folder, 'Fig3Groups_RhoUseOptZ_095.mat');
% resFileName = fullfile(folder, 'Fig3Groups_EVuseOptZ_095.mat');
resFileName = fullfile(folder, 'Fig3Groups_EVnoOptZ_095.mat');
% resFileName = fullfile(folder, 'comparisonRastersResMaxF.mat');
% resFileName = fullfile(folder, 'comparisonRastersRes.mat');
load(resFileName);
%% define some parameters
resFull = res;
resFull{end+1} = cell2mat(res);
zBinEdges = [0 110];
nZBins = numel(zBinEdges)-1;
rocThr = 0.95;
stdThr = 0;

%% plot rasters
axesLims = [-1 1];

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
        figName = sprintf('%s, z = [%1.0f %1.0f], roc < %4.2f', groupName{iGroup}, zBinEdges(iZBin), zBinEdges(iZBin+1), rocThr);
        figure('Name', figName, 'Position', [300 300 1200 400]);
        
%         idx = zPeak>=zBinEdges(iZBin) & zPeak < zBinEdges(iZBin+1);
%         idx = idx & rocArea<rocThr & thetaStd > stdThr;
        idx = ~isnan(rhoZTh);
        nSelectedCells = sum(idx);
        
        subplot(1, 3, 1)
        plot(rhoZTh(idx), rhoZD(idx), '.');
        hold on;
        plot(axesLims, axesLims, 'k:');
        axis equal
        xlim(axesLims)
        ylim(axesLims)
        plot([0 0], axesLims, 'k:', axesLims, [0 0], 'k:');
        xlabel('f(z, \theta)');
        ylabel('f(z, d)')
        box off;
        text(0.1, 0.9, sprintf('nCells = %1.0f', nSelectedCells), ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 16)

        
        subplot(1, 3, 2)
        plot(rhoZThD(idx), rhoZTh(idx), '.');
        hold on;
        plot(axesLims, axesLims, 'k:');
        axis equal
        xlim(axesLims)
        ylim(axesLims)
        plot([0 0], axesLims, 'k:', axesLims, [0 0], 'k:');
        xlabel('f(z, d, \theta)')
        ylabel('f(z, \theta)');
        title(figName);
        box off;
        
        
        subplot(1, 3, 3)
        plot(rhoZThD(idx), rhoZD(idx), '.');
        hold on;
        plot(axesLims, axesLims, 'k:');
        axis equal
        xlim(axesLims)
        ylim(axesLims)
        plot([0 0], axesLims, 'k:', axesLims, [0 0], 'k:');
        xlabel('f(z, d, \theta)');
        ylabel('f(z, d)')
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

