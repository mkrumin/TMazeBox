
folder = 'G:\DATA\';
% ExpRef = '2014-08-15_1931_MK012';
ExpRef = '2015-08-02_1709_MK020';
ExpRef = '2014-08-05_1937_MK014';
filename = fullfile(folder, [ExpRef, '_TMwExtras.mat']);
load(filename);
%%
EV = [];
planeIdx = [];
roiIdx = [];
roiClass = [];
for iPlane = TM.Planes
%     tmp = TM.modelEV(iPlane).ZThCV';
    tmp = [TM.modelFits{iPlane}.ZTh]';
    tmp = 1 - [tmp.errVal]';
    EV = [EV; tmp];
    nROIs = numel(tmp);
    planeIdx = cat(1, planeIdx, iPlane*ones(nROIs, 1));
    roiIdx = cat(1, roiIdx, [1:nROIs]');
    roiClass = cat(1, roiClass, cell2mat([TM.data2p{iPlane}.ROI.CellClasses]'));
end
planeIdx = planeIdx(roiClass == 's');
roiIdx = roiIdx(roiClass == 's');
EV = EV(roiClass == 's');

[~, sortedIdx] = sort(EV, 'descend');

%% session-specific analysis
nTrials = TM.dataTMaze.nTrials;
report = TM.dataTMaze.report;
idxR = report == 'R';
idxL = report == 'L';

nBins = 21;
for iPlane = TM.Planes
    [z, th] = binZTh({TM.modelExtras{iPlane}.z}', {TM.modelExtras{iPlane}.theta}', nBins);
    z = z(:);
    nTrials = size(th, 2);
    %             zz = repmat(z(:), 1, nTrials);
    % do some 'smart' interpolation to remove NaNs
    thFilt = th;
    for iTrial = 1:nTrials
        idxnan = isnan(th(:, iTrial));
        try
            thFilt(idxnan, iTrial) = interp1(z(~idxnan), th(~idxnan, iTrial), z(idxnan), 'linear');
            iFirst = find(~isnan(thFilt(:, iTrial)), 1, 'first');
            iLast = find(~isnan(thFilt(:, iTrial)), 1, 'last');
            thFilt(1:iFirst-1, iTrial) = thFilt(iFirst, iTrial);
            thFilt(iLast+1:end, iTrial) = thFilt(iLast, iTrial);
        catch
        end
    end
    th = thFilt;
    thStd{iPlane} = nanstd(th, [], 2);
    thROC{iPlane} = nan(nBins, 1);
    zROC{iPlane} = z;
    for iBin = 1:nBins
        thROC{iPlane}(iBin) = rocArea(th(iBin, idxL), th(iBin, idxR));
    end
end

nCells = length(EV);

%% cell-specific analysis

useModelForZOpt = false;

for iCell = 1:nCells
    fprintf('Cell %1.0f/%1.0f\n', iCell, nCells);
    iPlane = planeIdx(sortedIdx(iCell));
    iROI = roiIdx(sortedIdx(iCell));
    modelZTh = TM.modelFits{iPlane}(iROI).ZTh;
    modelZThD = TM.modelFits{iPlane}(iROI).ZThD;
    modelZD = TM.modelFits{iPlane}(iROI).ZD;
    
    zz = {TM.modelExtras{iPlane}.z}';
    th = {TM.modelExtras{iPlane}.theta}';
    oldZAxis = TM.trainingData{iPlane}(iROI).zThetaBinCentres{1};
    oldThAxis = TM.trainingData{iPlane}(iROI).zThetaBinCentres{2}; 
    oldMap = TM.trainingData{iPlane}(iROI).zThetaMap;
    zDMapR = TM.modelFits{iPlane}(iROI).ZD.zThetaMapR;
    zDMapL = TM.modelFits{iPlane}(iROI).ZD.zThetaMapL;
    
    thAxis = modelZTh.zThetaBinCentres{2};
    zAxis = modelZTh.zThetaBinCentres{1};
    
    
    for iTrial = 1:nTrials
        data{iTrial} = TM.modelExtras{iPlane}(iTrial).fData(:, iROI);
        
        fitZTh{iTrial} = TM.modelExtras{iPlane}(iTrial).fZThCV(:, iROI);
        fitZThD{iTrial} = TM.modelExtras{iPlane}(iTrial).fZThDCV(:, iROI);
        fitZD{iTrial} = TM.modelExtras{iPlane}(iTrial).fZDCV(:, iROI);
    end
    
    if useModelForZOpt
%         [~, zMaxInd] = max(abs(zDMapR-zDMapL));
        [~, zMaxInd] = max((zDMapR*sum(idxR) + zDMapL*sum(idxL))/nTrials);
        zPeak = zAxis(zMaxInd);
        zRange = zPeak + [-10 10];
        
    else
        [zUniform, dataUniform] = binZTh({TM.modelExtras{iPlane}.z}', data, 30);
        for iTrial = 1:nTrials
            idxnan = isnan(dataUniform(:, iTrial));
            try
                dataUniform(idxnan, iTrial) = interp1(zUniform(~idxnan), dataUniform(~idxnan, iTrial), zUniform(idxnan), 'linear');
            catch
            end
        end
        
        meanTrace = nanmean(dataUniform, 2);
        tooManyNans = sum(isnan(dataUniform), 2) >= nTrials/2;
        meanTrace(tooManyNans) = nan;
        try
            meanTrace = smooth(meanTrace, 'rlowess');
        end
        [~, zMaxInd] = max(meanTrace);
        zRange = zUniform(zMaxInd) + [-10 10];
        zPeak = zUniform(zMaxInd);
    end
    
    for iTrial = 1:nTrials
        idx = find(zz{iTrial}>=zRange(1) & zz{iTrial}<=zRange(2));
        fVector(iTrial) = nanmean(data{iTrial}(idx));
        thVector(iTrial) = nanmean(th{iTrial}(idx));
        zthVector(iTrial) = nanmean(fitZTh{iTrial}(idx));
        zthdVector(iTrial) = nanmean(fitZThD{iTrial}(idx));
        zdVector(iTrial) = nanmean(fitZD{iTrial}(idx));
    end
    
    idxValidTheta = abs(thVector) < max(abs(thAxis))-0.1;
    nThetaValidTrials = sum(idxValidTheta & ~isnan(fVector));
    nTotalTrials = sum(~isnan(fVector));
    zIdx = find(zAxis >= zRange(1) & zAxis <= zRange(2));
    zthFit = mean(modelZTh.zThetaMap(zIdx, :));
    zdFitR = mean(modelZD.zThetaMapR(zIdx))*ones(size(thAxis));
    zdFitL = mean(modelZD.zThetaMapL(zIdx))*ones(size(thAxis));
    zthdFitR = mean(modelZThD.zThetaMapR(zIdx, :));
    zthdFitL = mean(modelZThD.zThetaMapL(zIdx, :));
    
    rawMap = modelZTh.rawDataMap;
    rawMapR = modelZThD.rawDataMapR;
    rawMapL = modelZThD.rawDataMapL;
    
    zThMap = modelZTh.zThetaMap;
    zThMapR = modelZThD.zThetaMapR;
    zThMapL = modelZThD.zThetaMapL;
    zMapR = modelZD.zThetaMapR;
    zMapL = modelZD.zThetaMapL;
    
    clim = [min([rawMapR(:); rawMapL(:); rawMap(:)]), max([rawMapR(:); rawMapL(:); rawMap(:)])];
    
    nRows = 4;
    nColumns = 6;
    figName = sprintf('%s, iPlane = %1.0f, iROI = %1.0f', ExpRef, iPlane, iROI);
    hFig = figure('Name', figName, 'Position', [100 2 1600 1100]);
    colormap jet;
    
    subplot(nRows, nColumns, 1)
    imagesc(thAxis, zAxis, rawMapL);
    axis xy equal tight off
    caxis(clim);
    title('L trials');
    text(-100, 55, 'Data', 'FontSize', 20)
    hold on;
    plot(xlim, [zRange(1), zRange(1)], 'w--', 'LineWidth', 2)
    plot(xlim, [zRange(2), zRange(2)], 'w--', 'LineWidth', 2)
    
    subplot(nRows, nColumns, 2)
    imagesc(thAxis, zAxis, rawMapR);
    axis xy equal tight off
    caxis(clim);
    title('R trials');
    hold on;
    plot(xlim, [zRange(1), zRange(1)], 'w--', 'LineWidth', 2)
    plot(xlim, [zRange(2), zRange(2)], 'w--', 'LineWidth', 2)
    
    subplot(nRows, nColumns, 3)
    imagesc(thAxis, zAxis, rawMap);
    axis xy equal tight off
    caxis(clim);
    title('All trials');
    hold on;
    plot(xlim, [zRange(1), zRange(1)], 'w--', 'LineWidth', 2)
    plot(xlim, [zRange(2), zRange(2)], 'w--', 'LineWidth', 2)
    
    % plot trajectories
    subplot(nRows, nColumns, 4)
    for iTrial = 1:nTrials
        if idxR(iTrial)
            plot(th{iTrial}, zz{iTrial}, 'b');
        elseif idxL(iTrial)
            plot(th{iTrial}, zz{iTrial}, 'r');
        end
        hold on;
    end,
    axis equal
    xlim([-30 30]);
    ylim([5 105]);
    title('Trajectories')
    plot(xlim, [zRange(1), zRange(1)], 'k--', 'LineWidth', 2)
    plot(xlim, [zRange(2), zRange(2)], 'k--', 'LineWidth', 2)
    axis off
    
    
    subplot(nRows, nColumns, 5)
    plot(thROC{iPlane}, zROC{iPlane}, 'b', 'LineWidth', 3);
    hold on;
    %     axis equal
    xlim([-1 1]);
    ylim([5 105]);
    plot([0.9 0.9], ylim, ':');
    cellROC = interp1(zROC{iPlane}, thROC{iPlane}, zPeak);
    title(sprintf('AuROC_\\theta = %4.2f', cellROC));
    plot(xlim, [zRange(1), zRange(1)], 'k--', 'LineWidth', 2)
    plot(xlim, [zRange(2), zRange(2)], 'k--', 'LineWidth', 2)
    %     axis off
    
    subplot(nRows, nColumns, 7)
    imagesc(thAxis, zAxis, repmat(zMapL, 1, numel(thAxis)));
    axis xy equal tight off
    caxis(clim);
    text(-100, 55, 'f(z, d)', 'FontSize', 20)
    hold on;
    plot(xlim, [zRange(1), zRange(1)], 'w--', 'LineWidth', 2)
    plot(xlim, [zRange(2), zRange(2)], 'w--', 'LineWidth', 2)
    
    subplot(nRows, nColumns, 8)
    imagesc(thAxis, zAxis, repmat(zMapR, 1, numel(thAxis)));
    axis xy equal tight off
    caxis(clim);
    hold on;
    plot(xlim, [zRange(1), zRange(1)], 'w--', 'LineWidth', 2)
    plot(xlim, [zRange(2), zRange(2)], 'w--', 'LineWidth', 2)
    
    subplot(nRows, nColumns, 13)
    imagesc(thAxis, zAxis, zThMapL);
    axis xy equal tight off
    caxis(clim);
    text(-100, 55, 'f(z, \theta, d)', 'FontSize', 20)
    hold on;
    plot(xlim, [zRange(1), zRange(1)], 'w--', 'LineWidth', 2)
    plot(xlim, [zRange(2), zRange(2)], 'w--', 'LineWidth', 2)
    
    ax = subplot(nRows, nColumns, 14);
    imagesc(thAxis, zAxis, zThMapR);
    axis xy equal tight off
    caxis(clim);
    hold on;
    plot(xlim, [zRange(1), zRange(1)], 'w--', 'LineWidth', 2)
    plot(xlim, [zRange(2), zRange(2)], 'w--', 'LineWidth', 2)
    
    % this is a fake axis just to get the text in the right place
    subplot(nRows, nColumns, 19);
    % use limits of the previous axes
    xlim(ax.XLim);
    ylim(ax.YLim);
    axis xy equal off
    text(-100, 55, 'f(z, \theta)', 'FontSize', 20)
    
    subplot(nRows, nColumns, 21)
    imagesc(thAxis, zAxis, zThMap);
    axis xy equal tight off
    caxis(clim);
    hold on;
    plot(xlim, [zRange(1), zRange(1)], 'w--', 'LineWidth', 2)
    plot(xlim, [zRange(2), zRange(2)], 'w--', 'LineWidth', 2)
    
    supportR = thAxis>=prctile(thVector(idxR), 1) & thAxis<= prctile(thVector(idxR), 99);
    supportL = thAxis>=prctile(thVector(idxL), 1) & thAxis<= prctile(thVector(idxL), 99);
    
    thFakeValue = max(abs(thVector)) + 1.5;
    ax = subplot(nRows, nColumns, [10 11]);
    plot(thVector(idxR & idxValidTheta), fVector(idxR & idxValidTheta), '.b', 'MarkerSize', 12);
    hold on;
    plot(thVector(idxL & idxValidTheta), fVector(idxL & idxValidTheta), '.r', 'MarkerSize', 12);
    plot(sign(thVector(idxR & ~idxValidTheta))*thFakeValue, fVector(idxR & ~idxValidTheta), 'ob', 'MarkerSize', 4)
    plot(sign(thVector(idxL & ~idxValidTheta))*thFakeValue, fVector(idxL & ~idxValidTheta), 'or', 'MarkerSize', 4)
    %     plot(thVector(idxR), zdVector(idxR), 'ob', thVector(idxL), zdVector(idxL), 'or')
    plot(thAxis(supportR), zdFitR(supportR), 'b', thAxis(supportL), zdFitL(supportL), 'r', 'LineWidth', 2)
    xlim(round([min(thAxis), max(thAxis)]));
    %     xlabel('\theta [deg]');
    ylabel('\DeltaF/F');
    %     title(sprintf('z-d (E.V.=%5.3f)', TM.modelEV(iPlane).ZDCV(iROI)));
    title({'z-d'; sprintf('nAll = %1.0f, validTh = %1.0f', nTotalTrials, nThetaValidTrials)});
    ax.XTickLabel = '';
    ax.Clipping = 'off';
    
    ax = subplot(nRows, nColumns, [16 17]);
    plot(thVector(idxR & idxValidTheta), fVector(idxR & idxValidTheta), '.b', 'MarkerSize', 12);
    hold on;
    plot(thVector(idxL & idxValidTheta), fVector(idxL & idxValidTheta), '.r', 'MarkerSize', 12);
    plot(sign(thVector(idxR & ~idxValidTheta))*thFakeValue, fVector(idxR & ~idxValidTheta), 'ob', 'MarkerSize', 4)
    plot(sign(thVector(idxL & ~idxValidTheta))*thFakeValue, fVector(idxL & ~idxValidTheta), 'or', 'MarkerSize', 4)
    %     plot(thVector(idxR), zthdVector(idxR), 'ob', thVector(idxL), zthdVector(idxL), 'or')
    plot(thAxis(supportR), zthdFitR(supportR), 'b', thAxis(supportL), zthdFitL(supportL), 'r', 'LineWidth', 2)
    xlim(round([min(thAxis), max(thAxis)]))
    %     xlabel('\theta [deg]');
    ylabel('\DeltaF/F');
    %     title(sprintf('z-\\theta-d (E.V.=%5.3f)', TM.modelEV(iPlane).ZThDCV(iROI)));
    title('z-\theta-d');
    ax.XTickLabel = '';
    ax.Clipping = 'off';
    
    ax = subplot(nRows, nColumns, [22 23]);
    plot(thVector(idxR & idxValidTheta), fVector(idxR & idxValidTheta), '.b', 'MarkerSize', 12);
    hold on;
    plot(thVector(idxL & idxValidTheta), fVector(idxL & idxValidTheta), '.r', 'MarkerSize', 12);
    plot(sign(thVector(idxR & ~idxValidTheta))*thFakeValue, fVector(idxR & ~idxValidTheta), 'ob', 'MarkerSize', 4)
    plot(sign(thVector(idxL & ~idxValidTheta))*thFakeValue, fVector(idxL & ~idxValidTheta), 'or', 'MarkerSize', 4)
    %     plot(thVector, zthVector, 'ok')
    plot(thAxis, zthFit, 'k', 'LineWidth', 2)
    xlim(round([min(thAxis), max(thAxis)]))
    xlabel('\theta [deg]');
    ylabel('\DeltaF/F');
    %     title(sprintf('z-\\theta (E.V.=%5.3f)', TM.modelEV(iPlane).ZThCV(iROI)));
    title('z-\theta');
    ax.Clipping = 'off';
    
    
    %     figName = sprintf('Rasters %s, iPlane = %1.0f, iROI = %1.0f', ExpRef, iPlane, iROI);
    %     hFig2 = figure('Name', figName, 'Position', [200 100 300 900]);
    
    subplot(nRows, nColumns, 12);
    plot(fVector(idxR & idxValidTheta), zdVector(idxR & idxValidTheta), '.b', 'MarkerSize', 12);
    hold on;
    plot(fVector(idxR & ~idxValidTheta), zdVector(idxR & ~idxValidTheta), 'ob', 'MarkerSize', 4);
    plot(fVector(idxL & idxValidTheta), zdVector(idxL & idxValidTheta), '.r', 'MarkerSize', 12);
    plot(fVector(idxL & ~idxValidTheta), zdVector(idxL & ~idxValidTheta), 'or', 'MarkerSize', 4);
    axis equal
    xlim([0 max(fVector)]);
    ylim([0 max(fVector)]);
    plot(xlim, xlim, 'k:')
    %     xlabel('\DeltaF/F_{Data}');
    ylabel('\DeltaF/F_{Model}');
    tmp = corrcoef(fVector(idxValidTheta), zdVector(idxValidTheta), 'rows', 'complete');
    rho = tmp(2);
    tmp = corrcoef(fVector(:), zdVector(:), 'rows', 'complete');
    rho_all = tmp(2);
    title(sprintf('\\rho = %5.3f (\\rho_{all} %5.3f)', rho, rho_all));
    box off
    
    subplot(nRows, nColumns, 18);
    plot(fVector(idxR & idxValidTheta), zthdVector(idxR & idxValidTheta), '.b', 'MarkerSize', 12);
    hold on;
    plot(fVector(idxR & ~idxValidTheta), zthdVector(idxR & ~idxValidTheta), 'ob', 'MarkerSize', 4);
    plot(fVector(idxL & idxValidTheta), zthdVector(idxL & idxValidTheta), '.r', 'MarkerSize', 12);
    plot(fVector(idxL & ~idxValidTheta), zthdVector(idxL & ~idxValidTheta), 'or', 'MarkerSize', 4);
    axis equal
    xlim([0 max(fVector)]);
    ylim([0 max(fVector)]);
    plot(xlim, xlim, 'k:')
    %     xlabel('\DeltaF/F_{Data}');
    ylabel('\DeltaF/F_{Model}');
    tmp = corrcoef(fVector(idxValidTheta), zthdVector(idxValidTheta), 'rows', 'complete');
    rho = tmp(2);
    tmp = corrcoef(fVector(:), zthdVector(:), 'rows', 'complete');
    rho_all = tmp(2);
    title(sprintf('\\rho = %5.3f (\\rho_{all} %5.3f)', rho, rho_all));
    box off
    
    subplot(nRows, nColumns, 24);
    plot(fVector(idxR & idxValidTheta), zthVector(idxR & idxValidTheta), '.b', 'MarkerSize', 12);
    hold on;
    plot(fVector(idxR & ~idxValidTheta), zthVector(idxR & ~idxValidTheta), 'ob', 'MarkerSize', 4);
    plot(fVector(idxL & idxValidTheta), zthVector(idxL & idxValidTheta), '.r', 'MarkerSize', 12);
    plot(fVector(idxL & ~idxValidTheta), zthVector(idxL & ~idxValidTheta), 'or', 'MarkerSize', 4);
    axis equal
    xlim([0 max(fVector)]);
    ylim([0 max(fVector)]);
    plot(xlim, xlim, 'k:')
    xlabel('\DeltaF/F_{Data}');
    ylabel('\DeltaF/F_{Model}');
    tmp = corrcoef(fVector(idxValidTheta), zthVector(idxValidTheta), 'rows', 'complete');
    rho = tmp(2);
    tmp = corrcoef(fVector(:), zthVector(:), 'rows', 'complete');
    rho_all = tmp(2);
    title(sprintf('\\rho = %5.3f (\\rho_{all} %5.3f)', rho, rho_all));
    box off
    
        fileName = sprintf('%s_Cell#%03.0f(iPlane=%1.0f,iROI=%1.0f).png', ExpRef, iCell, iPlane, iROI);
        folderName = fullfile('G:\Reports\', ExpRef);
        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end
        print(hFig, fullfile(folderName, fileName), '-dpng')
    %
%         pause();
        close(hFig);
end