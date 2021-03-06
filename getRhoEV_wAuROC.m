function out = getRhoEV_wAuROC(TM)

rocThreshold = 0.95;

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
nTrials = TM.dataTMaze.nTrials;

%% session-specific AuROC curve estimation
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

%% session-specific analysis

rhoZTh = nan(nCells, 1);
rhoZThD = nan(nCells, 1);
rhoZD = nan(nCells, 1);
evZTh = nan(nCells, 1);
evZThD = nan(nCells, 1);
evZD = nan(nCells, 1);
evZTh_SbyS = nan(nCells, 1);
evZThD_SbyS = nan(nCells, 1);
evZD_SbyS = nan(nCells, 1);

for iCell = 1:nCells
    if mod(iCell,100)==1
        fprintf('Cell %1.0f/%1.0f\n', iCell, nCells);
    end
    iPlane = planeIdx(iCell);
    iROI = roiIdx(iCell);
    data = cell(nTrials, 1);
    fitZTh = cell(nTrials, 1);
    fitZThD = cell(nTrials, 1);
    fitZD = cell(nTrials, 1);
    
    for iTrial = 1:nTrials
        rocTrial = interp1(zROC{iPlane}, thROC{iPlane}, TM.modelExtras{iPlane}(iTrial).z);
        idxValid = rocTrial <= rocThreshold;
        data{iTrial} = TM.modelExtras{iPlane}(iTrial).fData(idxValid, iROI);
        fitZTh{iTrial} = TM.modelExtras{iPlane}(iTrial).fZThCV(idxValid, iROI);
        fitZThD{iTrial} = TM.modelExtras{iPlane}(iTrial).fZThDCV(idxValid, iROI);
        fitZD{iTrial} = TM.modelExtras{iPlane}(iTrial).fZDCV(idxValid, iROI);
    end
    
    fVector = nan(nTrials, 1);
    zthVector = nan(nTrials, 1);
    zthdVector = nan(nTrials, 1);
    zdVector = nan(nTrials, 1);
    
    for iTrial = 1:nTrials
        fVector(iTrial) = nanmean(data{iTrial});
        zthVector(iTrial) = nanmean(fitZTh{iTrial});
        zthdVector(iTrial) = nanmean(fitZThD{iTrial});
        zdVector(iTrial) = nanmean(fitZD{iTrial});
    end
    
        tmp = corrcoef(fVector(:), zdVector(:), 'rows', 'complete');
        rhoZD(iCell) = tmp(2);
        
        tmp = corrcoef(fVector(:), zthdVector(:), 'rows', 'complete');
        rhoZThD(iCell) = tmp(2);
        
        tmp = corrcoef(fVector(:), zthVector(:), 'rows', 'complete');
        rhoZTh(iCell) = tmp(2);
        
        evZD(iCell) = 1 - nanmean((fVector(:)-zdVector(:)).^2)/...
            nanmean((fVector(:)-nanmean(fVector(:))).^2);
        
        evZThD(iCell) = 1 - nanmean((fVector(:)-zthdVector(:)).^2)/...
            nanmean((fVector(:)-nanmean(fVector(:))).^2);
        
        evZTh(iCell) = 1 - nanmean((fVector(:)-zthVector(:)).^2)/...
            nanmean((fVector(:)-nanmean(fVector(:))).^2);
        
        dataAll = cell2mat(data);
        fitZThAll = cell2mat(fitZTh);
        fitZThDAll = cell2mat(fitZThD);
        fitZDAll = cell2mat(fitZD);

        evZD_SbyS(iCell) = 1 - nanmean((dataAll(:)-fitZDAll(:)).^2)/...
            nanmean((dataAll(:)-nanmean(dataAll(:))).^2);
        
        evZThD_SbyS(iCell) = 1 - nanmean((dataAll(:)-fitZThDAll(:)).^2)/...
            nanmean((dataAll(:)-nanmean(dataAll(:))).^2);
        
        evZTh_SbyS(iCell) = 1 - nanmean((dataAll(:)-fitZThAll(:)).^2)/...
            nanmean((dataAll(:)-nanmean(dataAll(:))).^2);

end

out.ExpRef = TM.expRef;
out.rhoZTh = rhoZTh;
out.rhoZThD = rhoZThD;
out.rhoZD = rhoZD;
out.evZTh_TbyT = evZTh;
out.evZThD_TbyT = evZThD;
out.evZD_TbyT = evZD;
out.evZTh_SbyS = evZTh_SbyS;
out.evZThD_SbyS = evZThD_SbyS;
out.evZD_SbyS = evZD_SbyS;
out.planeIdx = planeIdx;
out.roiIdx = roiIdx;
% out.EV = EV;
