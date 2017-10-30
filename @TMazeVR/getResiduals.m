function getResiduals(obj)

for iPlane = obj.Planes
    fData = obj.data2p{iPlane}.F;
    tData = obj.times2p{iPlane};
    f0 = prctile(fData, 20);
    fData = bsxfun(@rdivide, bsxfun(@minus, fData, f0), f0);
    
    nCells = size(fData, 2);
    nTrials = obj.dataTMaze.nTrials;
    options.econ = true;
    for iTrial = 1:nTrials
        [thVector, zVector, fVector, ~] = buildVectors(obj, iTrial, tData, fData, options);
        residuals = nan(size(fVector));
        if ~isempty(thVector)
            for iCell = 1:nCells
                zAxis = obj.trainingData{iPlane}(iCell).zThetaBinCentres{1};
                thAxis = obj.trainingData{iPlane}(iCell).zThetaBinCentres{2};
%                 thAxis(1) = thAxis(1)-0.001;
%                 thAxis(end) = thAxis(end)+0.001;
                theMap = obj.trainingData{iPlane}(iCell).zThetaMap;
                fModel = interp2(thAxis, zAxis, theMap, thVector, zVector, 'spline');%, 0);
                residuals(:, iCell) = fVector(:, iCell) - fModel;
            end
        end
        obj.residualData{iPlane}(iTrial).residuals = residuals;
        obj.residualData{iPlane}(iTrial).th = thVector;
        obj.residualData{iPlane}(iTrial).z = zVector;
    end
    
    thVector = cell2mat({obj.residualData{iPlane}(:).th}');
    zVector = cell2mat({obj.residualData{iPlane}(:).z}');
    resVector = cell2mat({obj.residualData{iPlane}(:).residuals}');
    
    for iCell = 1:nCells
        zEdges = obj.trainingData{iPlane}(iCell).zEdges;
        thEdges = obj.trainingData{iPlane}(iCell).thEdges;
        optStd = obj.trainingData{iPlane}(iCell).optStd;
        hFilter = ndGaussian(optStd);
        [resMap, ] = estMap([zVector, thVector], resVector(:, iCell), {zEdges, thEdges}, hFilter);
        obj.trainingData{iPlane}(iCell).residualMap = resMap;
    end

end

% make it a column cell array
obj.residualData = obj.residualData(:);
