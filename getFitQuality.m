function out = getFitQuality(TM)

nBinsROC = 21;
nZBinsF = 30;
useModelForZOpt = false;
useZOpt = false;
rocThresh = 1;
useRho = false; % if false, use explained variance

EV = [];
planeIdx = [];
roiIdx = [];
roiClass = [];
for iPlane = TM.Planes
    %     tmp = TM.modelEV(iPlane).ZThCV';
    tmp = [TM.modelFits{iPlane}.ZTh]';
    nROIs = TM.nROIs(iPlane);
%     EV = [EV; tmp];
    planeIdx = cat(1, planeIdx, iPlane*ones(nROIs, 1));
    roiIdx = cat(1, roiIdx, [1:nROIs]');
    roiClass = cat(1, roiClass, cell2mat([TM.data2p{iPlane}.ROI.CellClasses]'));
end
idxS = roiClass == 's';
planeIdx = planeIdx(idxS);
roiIdx = roiIdx(idxS);
% EV = EV(idxS);

nCells = sum(idxS);

%% session-specific analysis
nTrials = TM.dataTMaze.nTrials;
idxR = TM.dataTMaze.report == 'R';
idxL = TM.dataTMaze.report == 'L';

for iPlane = TM.Planes
    [z, th] = binZTh({TM.modelExtras{iPlane}.z}', {TM.modelExtras{iPlane}.theta}', nBinsROC);
    z = z(:);
    %     nTrials = size(th, 2);
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
            % will get here if the whole th vector is NaNs
        end
    end
    th = thFilt;
    thStd{iPlane} = nanstd(th, [], 2);
    thROC{iPlane} = nan(nBinsROC, 1);
    zROC{iPlane} = z;
    for iBin = 1:nBinsROC
        thROC{iPlane}(iBin) = rocArea(th(iBin, idxL), th(iBin, idxR));
    end
end

%% cell-specific analysis

rhoZTh = nan(nCells, 1);
rhoZThD = nan(nCells, 1);
rhoZD = nan(nCells, 1);
rhoZTh_all = nan(nCells, 1);
rhoZThD_all = nan(nCells, 1);
rhoZD_all = nan(nCells, 1);
zPeak = nan(nCells, 1);

for iCell = 1:nCells
    if mod(iCell,100)==1
        fprintf('Cell %1.0f/%1.0f\n', iCell, nCells);
    end
    iPlane = planeIdx(iCell);
    iROI = roiIdx(iCell);
    modelZTh = TM.modelFits{iPlane}(iROI).ZTh;
    %     modelZThD = TM.modelFits{iPlane}(iROI).ZThD;
    %     modelZD = TM.modelFits{iPlane}(iROI).ZD;
    
    zz = {TM.modelExtras{iPlane}.z}';
    th = {TM.modelExtras{iPlane}.theta}';
    %     oldZAxis = TM.trainingData{iPlane}(iROI).zThetaBinCentres{1};
    %     oldThAxis = TM.trainingData{iPlane}(iROI).zThetaBinCentres{2};
    %     oldMap = TM.trainingData{iPlane}(iROI).zThetaMap;
    
    thAxis = modelZTh.zThetaBinCentres{2};
    zAxis = modelZTh.zThetaBinCentres{1};
    
    
    for iTrial = 1:nTrials
        data{iTrial} = TM.modelExtras{iPlane}(iTrial).fData(:, iROI);
        
        fitZTh{iTrial} = TM.modelExtras{iPlane}(iTrial).fZThCV(:, iROI);
        fitZThD{iTrial} = TM.modelExtras{iPlane}(iTrial).fZThDCV(:, iROI);
        fitZD{iTrial} = TM.modelExtras{iPlane}(iTrial).fZDCV(:, iROI);
    end
    
    if useModelForZOpt
        zDMapR = TM.modelFits{iPlane}(iROI).ZD.zThetaMapR;
        zDMapL = TM.modelFits{iPlane}(iROI).ZD.zThetaMapL;
        
        %         [~, zMaxInd] = max(abs(zDMapR-zDMapL));
        [~, zMaxInd] = max((zDMapR*sum(idxR) + zDMapL*sum(idxL))/nTrials);
        zPeak(iCell) = zAxis(zMaxInd);
    else
        [zUniform, dataUniform] = binZTh({TM.modelExtras{iPlane}.z}', data, nZBinsF);
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
        zPeak(iCell) = zUniform(zMaxInd);
    end
    
    zRange = zPeak + [-10 10];
    
    for iTrial = 1:nTrials
        if useZOpt
            idxInRange = zz{iTrial}>=zRange(1) & zz{iTrial}<=zRange(2);
        else
            idxInRange = true(size(zz{iTrial}));
        end    
        rocTmp = interp1(zROC{iPlane}, thROC{iPlane}, zz{iTrial});
        idxValidROC = rocTmp <= rocThresh;
        idx = idxInRange & idxValidROC;
        fVector(iTrial) = nanmean(data{iTrial}(idx));
        thVector(iTrial) = nanmean(th{iTrial}(idx));
        zthVector(iTrial) = nanmean(fitZTh{iTrial}(idx));
        zthdVector(iTrial) = nanmean(fitZThD{iTrial}(idx));
        zdVector(iTrial) = nanmean(fitZD{iTrial}(idx));
    end
    
    idxValidTheta = abs(thVector) < max(abs(thAxis))-0.1;
    idxValidTrials = idxValidTheta & ~isnan(fVector);
    nThetaValidTrials = sum(idxValidTrials);
    nTotalTrials = sum(~isnan(fVector));
    
    if useRho
        tmp = corrcoef(fVector(idxValidTrials), zdVector(idxValidTrials), 'rows', 'complete');
        rhoZD(iCell) = tmp(2);
        tmp = corrcoef(fVector(:), zdVector(:), 'rows', 'complete');
        rhoZD_all(iCell) = tmp(2);
        
        tmp = corrcoef(fVector(idxValidTrials), zthdVector(idxValidTrials), 'rows', 'complete');
        rhoZThD(iCell) = tmp(2);
        tmp = corrcoef(fVector(:), zthdVector(:), 'rows', 'complete');
        rhoZThD_all(iCell) = tmp(2);
        
        tmp = corrcoef(fVector(idxValidTrials), zthVector(idxValidTrials), 'rows', 'complete');
        rhoZTh(iCell) = tmp(2);
        tmp = corrcoef(fVector(:), zthVector(:), 'rows', 'complete');
        rhoZTh_all(iCell) = tmp(2);
    else
        
        rhoZD(iCell) = 1 - nanmean((fVector(idxValidTrials)-zdVector(idxValidTrials)).^2)/...
            nanmean((fVector(idxValidTrials)-nanmean(fVector(idxValidTrials))).^2);
        rhoZD_all(iCell) = 1 - nanmean((fVector(:)-zdVector(:)).^2)/...
            nanmean((fVector(:)-nanmean(fVector(:))).^2);
        
        rhoZThD(iCell) = 1 - nanmean((fVector(idxValidTrials)-zthdVector(idxValidTrials)).^2)/...
            nanmean((fVector(idxValidTrials)-nanmean(fVector(idxValidTrials))).^2);
        rhoZThD_all(iCell) = 1 - nanmean((fVector(:)-zthdVector(:)).^2)/...
            nanmean((fVector(:)-nanmean(fVector(:))).^2);
        
        rhoZTh(iCell) = 1 - nanmean((fVector(idxValidTrials)-zthVector(idxValidTrials)).^2)/...
            nanmean((fVector(idxValidTrials)-nanmean(fVector(idxValidTrials))).^2);
        rhoZTh_all(iCell) = 1 - nanmean((fVector(:)-zthVector(:)).^2)/...
            nanmean((fVector(:)-nanmean(fVector(:))).^2);
        
    end
    
    %     zIdx = find(zAxis >= zRange(1) & zAxis <= zRange(2));
    %     zthFit = mean(modelZTh.zThetaMap(zIdx, :));
    %     zdFitR = mean(modelZD.zThetaMapR(zIdx))*ones(size(thAxis));
    %     zdFitL = mean(modelZD.zThetaMapL(zIdx))*ones(size(thAxis));
    %     zthdFitR = mean(modelZThD.zThetaMapR(zIdx, :));
    %     zthdFitL = mean(modelZThD.zThetaMapL(zIdx, :));
    %
    %     rawMap = modelZTh.rawDataMap;
    %     rawMapR = modelZThD.rawDataMapR;
    %     rawMapL = modelZThD.rawDataMapL;
    %
    %     zThMap = modelZTh.zThetaMap;
    %     zThMapR = modelZThD.zThetaMapR;
    %     zThMapL = modelZThD.zThetaMapL;
    %     zMapR = modelZD.zThetaMapR;
    %     zMapL = modelZD.zThetaMapL;
    
end

out.ExpRef = TM.expRef;
out.rhoZTh = rhoZTh;
out.rhoZThD = rhoZThD;
out.rhoZD = rhoZD;
out.rhoZTh_all = rhoZTh_all;
out.rhoZThD_all = rhoZThD_all;
out.rhoZD_all = rhoZD_all;
out.zPeak = zPeak;
out.thStd = thStd;
out.thROC = thROC;
out.zROC = zROC;
out.planeIdx = planeIdx;
out.roiIdx = roiIdx;
% out.EV = EV;
