function trainSFZDModels(obj, iPlane, iROI, options)

doPlotting = false; % flag for debugging
if nargin<4
    options = struct;
end
options = fillOptions(obj, options);

fData = obj.data2p{iPlane}.F(:, iROI);

F0 = prctile(fData, 20);
fData = (fData-F0)/F0;

tData = obj.times2p{iPlane};
tData = tData - options.delay;

trialIdx = 1:length(obj.dataTMaze.contrastSequence);
nTrials = length(trialIdx);
thAll = cell(nTrials, 1);
zAll = cell(nTrials, 1);
fAll = cell(nTrials, 1);
validTrials = true(nTrials, 1);
for iTrial = 1:nTrials
    [thAll{iTrial}, zAll{iTrial}, fAll{iTrial}] = buildVectors(obj, iTrial, tData, fData, options);
    validTrials(iTrial) = ~isempty(fAll{iTrial});
end

zMin = min(abs(cell2mat(zAll)));
zMax = max(abs(cell2mat(zAll)));
thMax = max(abs(cell2mat(thAll)));
% clear zAll thAll;

% now, a complicated way to define bin edges
thEdges = buildEdges(thMax, options.dTheta);
zEdges = buildEdges([zMin, zMax], options.dZ);

nChunks = options.cvFactor;
validTrials = validTrials & ismember(obj.dataTMaze.report, 'LR')';
validTrialIdx = find(validTrials);
report = obj.dataTMaze.report(validTrials);
report = report(:);

cvObject = cvpartition(report, 'kFold', nChunks);
fValid = fAll(validTrials);
zValid = zAll(validTrials);
thValid = thAll(validTrials);

thVector = cell(nChunks, 1);
zVector = cell(nChunks, 1);
fVector = cell(nChunks, 1);
zR = cell(nChunks, 1);
thR = cell(nChunks, 1);
fR = cell(nChunks, 1);
zL = cell(nChunks, 1);
thL = cell(nChunks, 1);
fL = cell(nChunks, 1);
coords = cell(nChunks, 1);
coordsR = cell(nChunks, 1);
coordsL = cell(nChunks, 1);
for iChunk = 1:nChunks
    idxL = cvObject.test(iChunk) & report == 'L';
    idxR = cvObject.test(iChunk) & report == 'R';
    thVector{iChunk} = cell2mat(thValid([idxL | idxR]));
    zVector{iChunk} = cell2mat(zValid([idxL | idxR]));
    fVector{iChunk} = cell2mat(fValid([idxL | idxR]));
    zR{iChunk} = cell2mat(zValid(idxR));
    thR{iChunk} = cell2mat(thValid(idxR));
    fR{iChunk} = cell2mat(fValid(idxR));
    zL{iChunk} = cell2mat(zValid(idxL));
    thL{iChunk} = cell2mat(thValid(idxL));
    fL{iChunk} = cell2mat(fValid(idxL));
    coords{iChunk} = [zVector{iChunk}, thVector{iChunk}];
    coordsR{iChunk} = [zR{iChunk}, thR{iChunk}];
    coordsL{iChunk} = [zL{iChunk}, thL{iChunk}];
end

% estimating the optimal filter STDs
[optStd, errVals, errValMatrix, stdGridValues] = ...
    estOptimalStd(coords, fVector, {zEdges, thEdges}, options.errPower);
[optStdR3D, errValR3D, errValMatrixR3D, stdGridValuesR3D] = ...
    estOptimalStd(coordsR, fR, {zEdges, thEdges}, options.errPower);
[optStdL3D, errValL3D, errValMatrixL3D, stdGridValuesL3D] = ...
    estOptimalStd(coordsL, fL, {zEdges, thEdges}, options.errPower);
[optStdR, errValR, errValMatrixR, stdGridValuesR] = ...
    estOptimalStd1D(zR, fR, {zEdges}, options.errPower);
[optStdL, errValL, errValMatrixL, stdGridValuesL] = ...
    estOptimalStd1D(zL, fL, {zEdges}, options.errPower);

% [optStd25D, errVals25D, ~, ~] = ...
%     estOptimalStd25D([coordsL, coordsR], [fL, fR], {zEdges, thEdges}, options.errPower);
% [optStd15D, errVals15D, ~, ~] = ...
%     estOptimalStd15D([zL, zR], [fL, fR], {zEdges}, options.errPower);

nR = sum(report == 'R');
nL = sum(report == 'L');

% these are the results of the main z-theta model
obj.modelFits{iPlane}(iROI).ZTh.optStd = optStd;
obj.modelFits{iPlane}(iROI).ZTh.errVal = errVals;
obj.modelFits{iPlane}(iROI).ZTh.errValMatrix = errValMatrix;
obj.modelFits{iPlane}(iROI).ZTh.stdGridValues = stdGridValues;
obj.modelFits{iPlane}(iROI).ZTh.zEdges = zEdges;
obj.modelFits{iPlane}(iROI).ZTh.thEdges = thEdges;

hFilter = ndGaussian(optStd);
[theMap, binCentres] = estMap(cell2mat(coords), cell2mat(fVector), {zEdges, thEdges}, hFilter);
hFilter = ndGaussian([0.2 0.2]);
[rawMap, ~] = estMap(cell2mat(coords), cell2mat(fVector), {zEdges, thEdges}, hFilter);

obj.modelFits{iPlane}(iROI).ZTh.options = options;
obj.modelFits{iPlane}(iROI).ZTh.zThetaMap = theMap;
obj.modelFits{iPlane}(iROI).ZTh.zThetaBinCentres = binCentres;
obj.modelFits{iPlane}(iROI).ZTh.rawDataMap = rawMap;

% these are the results of the full z-theta-d model
% (z-theta maps for each d can be completely different)
obj.modelFits{iPlane}(iROI).ZThD.optStdR = optStdR3D;
obj.modelFits{iPlane}(iROI).ZThD.optStdL = optStdL3D;
obj.modelFits{iPlane}(iROI).ZThD.errValR = errValR3D;
obj.modelFits{iPlane}(iROI).ZThD.errValL = errValL3D;
obj.modelFits{iPlane}(iROI).ZThD.errVal = (errValL3D*nL + errValR3D*nR)/(nL + nR);
obj.modelFits{iPlane}(iROI).ZThD.errValMatrixR = errValMatrixR3D;
obj.modelFits{iPlane}(iROI).ZThD.errValMatrixL = errValMatrixL3D;
obj.modelFits{iPlane}(iROI).ZThD.stdGridValuesR = stdGridValuesR3D;
obj.modelFits{iPlane}(iROI).ZThD.stdGridValuesL = stdGridValuesL3D;
obj.modelFits{iPlane}(iROI).ZThD.zEdges = zEdges;
obj.modelFits{iPlane}(iROI).ZThD.thEdges = thEdges;
obj.modelFits{iPlane}(iROI).ZThD.options = options;

hFilter = ndGaussian(optStdR3D);
[theMapR, binCentres] = estMap(cell2mat(coordsR), cell2mat(fR), {zEdges, thEdges}, hFilter);
hFilter = ndGaussian(optStdL3D);
[theMapL, ~] = estMap(cell2mat(coordsL), cell2mat(fL), {zEdges, thEdges}, hFilter);
hFilter = ndGaussian([0.2 0.2]);
[rawMapR, ~] = estMap(cell2mat(coordsR), cell2mat(fR), {zEdges, thEdges}, hFilter);
[rawMapL, ~] = estMap(cell2mat(coordsL), cell2mat(fL), {zEdges, thEdges}, hFilter);

obj.modelFits{iPlane}(iROI).ZThD.zThetaMapR = theMapR;
obj.modelFits{iPlane}(iROI).ZThD.zThetaMapL = theMapL;
obj.modelFits{iPlane}(iROI).ZThD.zThetaBinCentres = binCentres;
obj.modelFits{iPlane}(iROI).ZThD.rawDataMapR = rawMapR;
obj.modelFits{iPlane}(iROI).ZThD.rawDataMapL = rawMapL;

%%

if doPlotting
    figure
    clim = [min([theMapR(:); theMapL(:); theMap(:)]), max([theMapR(:); theMapL(:); theMap(:)])];
    clim = [min([rawMapR(:); rawMapL(:); rawMap(:)]), max([rawMapR(:); rawMapL(:); rawMap(:)])];
    subplot(2, 3, 1);
    imagesc(binCentres{2}, binCentres{1}, rawMapL);
    axis xy equal tight;
    caxis(clim);
    subplot(2, 3, 2);
    imagesc(binCentres{2}, binCentres{1}, rawMap);
    axis xy equal tight;
    caxis(clim);
    subplot(2, 3, 3);
    imagesc(binCentres{2}, binCentres{1}, rawMapR);
    axis xy equal tight;
    caxis(clim);
    
    subplot(2, 3, 4);
    imagesc(binCentres{2}, binCentres{1}, theMapL);
    axis xy equal tight;
    caxis(clim);
    subplot(2, 3, 5);
    imagesc(binCentres{2}, binCentres{1}, theMap);
    axis xy equal tight;
    caxis(clim);
    subplot(2, 3, 6);
    imagesc(binCentres{2}, binCentres{1}, theMapR);
    axis xy equal tight;
    caxis(clim);
end

%%
% these are the results of the 'full' z-d model
% (z tuning for each d can be completely different)
obj.modelFits{iPlane}(iROI).ZD.optStdR = optStdR;
obj.modelFits{iPlane}(iROI).ZD.optStdL = optStdL;
obj.modelFits{iPlane}(iROI).ZD.errValR = errValR;
obj.modelFits{iPlane}(iROI).ZD.errValL = errValL;
obj.modelFits{iPlane}(iROI).ZD.errVal = (errValL*nL + errValR*nR)/(nL + nR);
obj.modelFits{iPlane}(iROI).ZD.errValMatrixR = errValMatrixR;
obj.modelFits{iPlane}(iROI).ZD.errValMatrixL = errValMatrixL;
obj.modelFits{iPlane}(iROI).ZD.stdGridValuesR = stdGridValuesR;
obj.modelFits{iPlane}(iROI).ZD.stdGridValuesL = stdGridValuesL;
obj.modelFits{iPlane}(iROI).ZD.zEdges = zEdges;
obj.modelFits{iPlane}(iROI).ZD.options = options;


hFilter = ndGaussian(optStdR);
[theMapR, binCentres] = estMap(cell2mat(zR), cell2mat(fR), {zEdges}, hFilter);
hFilter = ndGaussian(optStdL);
[theMapL, ~] = estMap(cell2mat(zL), cell2mat(fL), {zEdges}, hFilter);
hFilter = ndGaussian(0.2);
[rawMapR, ~] = estMap(cell2mat(zR), cell2mat(fR), {zEdges}, hFilter);
[rawMapL, ~] = estMap(cell2mat(zL), cell2mat(fL), {zEdges}, hFilter);

obj.modelFits{iPlane}(iROI).ZD.zThetaMapR = theMapR;
obj.modelFits{iPlane}(iROI).ZD.zThetaMapL = theMapL;
obj.modelFits{iPlane}(iROI).ZD.zThetaBinCentres = binCentres;
obj.modelFits{iPlane}(iROI).ZD.rawDataMapR = rawMapR;
obj.modelFits{iPlane}(iROI).ZD.rawDataMapL = rawMapL;

%%
if doPlotting
    figure
    zAxis = binCentres{1};
    plot(zAxis, theMapR, zAxis, theMapL, 'LineWidth', 2)
    hold on;
    plot(zAxis, rawMapR, '--', zAxis, rawMapL, '--', 'LineWidth', 1)
    xlim([min(zEdges), max(zEdges)])
    legend('Model_R', 'Model_L', 'Raw_R', 'Raw_L');
    xlabel('z [cm]');
    ylabel('\DeltaF/F');
end
end % TrainMaps()

%==========================================================================
function optOut = fillOptions(obj, optIn)

optOut = optIn;

if ~isfield(optIn, 'dZ')
    optOut.dZ = 3; % [cm]
end

if ~isfield(optIn, 'dTheta')
    optOut.dTheta = 2; % [deg]
end

if ~isfield(optIn, 'cvFactor')
    optOut.cvFactor = 10; %
end

if ~isfield(optIn, 'delay')
    optOut.delay = 0; %
end

if ~isfield(optIn, 'errPower')
    optOut.errPower = 2; %
end

if ~isfield(optIn, 'econ')
    optOut.econ = true; %
end

end %fillOptions();

%==========================================================================
function [stdOut, errVal, errValMatrix, gridValues] = estOptimalStd(coords, signal, binEdges, errPow)

coords = coords(:);
signal = signal(:);
nTrials = length(coords);
nDims = size(coords{1}, 2);

% precalculate the occupancy and signal maps (unsmoothed - filter will be
% our optimization parameter),
% also linear indices of the testCoordinates
nRows = length(binEdges{1})-1;
nColumns = length(binEdges{2})-1;
occMaps = zeros(nRows, nColumns, nTrials);
fMaps = zeros(nRows, nColumns, nTrials);
testXYLinearIdx = cell(nTrials, 1);
for iTrial = 1:nTrials
    if isempty(coords{iTrial})
        continue;
    end
    
    [occMaps(:,:,iTrial), binCentres] = buildOccupMap(coords{iTrial}, binEdges);
    fMaps(:,:,iTrial) = buildAccumMap(coords{iTrial}, signal{iTrial}, binEdges);
    ind = nan(size(coords{iTrial}));
    [~, ind(:, 1)] = histc(coords{iTrial}(:,1), binEdges{1});
    [~, ind(:, 2)] = histc(coords{iTrial}(:,2), binEdges{2});
    testXYLinearIdx{iTrial} = sub2ind(size(occMaps(:,:,1)), ind(:,1), ind(:,2));
end

occMapsXOR = zeros(size(occMaps));
fMapsXOR = zeros(size(fMaps));
for iTrial = 1:nTrials
    idx = setdiff(1:nTrials, iTrial);
    occMapsXOR(:,:,iTrial) = sum(occMaps(:,:,idx), 3);
    fMapsXOR(:,:,iTrial) = sum(fMaps(:,:,idx), 3);
end

% let's run a grid search for the optimal parameters

sz = size(occMaps);
gridValues = cell(nDims, 1);
for iDim = 1:nDims
    minValue = 1;
    maxValue = sz(iDim);
    pp = 3;
    gridValues{iDim} = linspace(minValue^(1/pp), maxValue^(1/pp), 15).^pp;
end

if nDims == 2
    errValMatrix = nan(length(gridValues{1}), length(gridValues{2}));
    for nInd = 1:length(gridValues{2})
        % prefilter along the 2nd dimension
        x2 = gridValues{2}(nInd);
        hGauss2 = ndGaussian(x2);
        
        meanVal = squeeze(sum(sum(fMapsXOR), 2)./sum(sum(occMapsXOR), 2));
        
        occMapsXORSpecial = permute(occMapsXOR, [2, 1, 3]);
        occMapsXORSpecial = reshape(occMapsXORSpecial, nColumns, nRows*nTrials);
        occMapsXORFilt = conv2(hGauss2, 1, occMapsXORSpecial, 'same');
        occMapsXORFilt = reshape(occMapsXORFilt, nColumns, nRows, nTrials);
        occMapsXORFilt = permute(occMapsXORFilt, [2, 1, 3]);
        
        fMapsXORSpecial = permute(fMapsXOR, [2, 1, 3]);
        fMapsXORSpecial = reshape(fMapsXORSpecial, nColumns, nRows*nTrials);
        fMapsXORFilt = conv2(hGauss2, 1, fMapsXORSpecial, 'same');
        fMapsXORFilt = reshape(fMapsXORFilt, nColumns, nRows, nTrials);
        fMapsXORFilt = permute(fMapsXORFilt, [2, 1, 3]);
        
        for mInd = 1:length(gridValues{1})
            x1 = gridValues{1}(mInd);
            errValMatrix(mInd, nInd) = mapError(x1, occMapsXORFilt, fMapsXORFilt, meanVal, testXYLinearIdx, signal, errPow);
        end
    end
else
    fprint('For now, grid search is implemented only for 2D grids\n');
end

ind = find(errValMatrix == min(errValMatrix(:)));
[mOpt, nOpt] = ind2sub(size(errValMatrix), ind);
stdOut = [gridValues{1}(mOpt), gridValues{2}(nOpt)];
errVal = errValMatrix(ind);

% % now run the optimization
% initStd = 5*ones(1, nDims);
% [stdOut, errVal] = fminsearch(@(x) mapError(x, occMaps, fMaps, testXYLlinIdx, testSignal), initStd);

end % estOptimalStd()

%==========================================================================
function errValue = mapError(x1, occMaps, fMaps, meanF, testXYIdx, testF, errPow)

x1(x1<0.2) = 0.2; % HACK! This value is in pixels, not real units

hGauss = ndGaussian(x1);

[nRows, nColumns, nTrials] = size(occMaps);
chunkErrors = nan(nTrials, 1);

epsilon = 0.01;
meanFMatrix = repmat(permute(meanF, [2, 3, 1]), nRows, nColumns, 1);
% allPredictionMaps = (imfilter(fMaps, hGauss)+epsilon*meanFMatrix)./(imfilter(occMaps, hGauss)+epsilon);
fMapsSpecial = reshape(fMaps, nRows, []);
occMapsSpecial = reshape(occMaps, nRows, []);
meanFSpecial = reshape(meanFMatrix, nRows, []);
allPredictionMaps = (conv2(hGauss, 1, fMapsSpecial, 'same')+epsilon*meanFSpecial)./(conv2(hGauss, 1, occMapsSpecial, 'same')+epsilon);
allPredictionMaps = reshape(allPredictionMaps, nRows, nColumns, nTrials);

for iTrial = 1:nTrials
    if isempty(testF{iTrial})
        % there is no data for this trial
        continue;
    end
    %     occM = occMaps(:,:,iTrial);
    %     fM = fMaps(:,:,iTrial);
    %     predictionMap = filterAndDivideMaps1D(occM, fM, meanF(iTrial), hGauss);
    predictionMap = allPredictionMaps(:,:,iTrial);
    
    predictedF = predictionMap(testXYIdx{iTrial});
    
    %     chunkErrors(iChunk) = var(predictedF - testF{iChunk})/var(testF{iChunk});
    meanFTraining = meanF(iTrial);
    % the idea here is to check how much better the map is relative to just
    % using mean F map.
    switch errPow
        case 1
            chunkErrors(iTrial) = sum(abs(predictedF - testF{iTrial}))/sum(abs(meanFTraining - testF{iTrial}));
        case 2
            chunkErrors(iTrial) = sum((predictedF - testF{iTrial}).^2)/sum((meanFTraining - testF{iTrial}).^2);
        otherwise
            chunkErrors(iTrial) = sum(abs(predictedF - testF{iTrial}).^errPow)/sum(abs(meanFTraining - testF{iTrial}).^errPow);
    end
end
% errValue = nanmedian(chunkErrors);
errValue = nanmean(chunkErrors); % mean() works better, maps are smoother and nicer, "outliers have less effect"

end % mapErrorFaster()

%==========================================================================
function [stdOut, errVal, errValMatrix, gridValues] = estOptimalStd1D(coords, signal, binEdges, errPow)

coords = coords(:);
signal = signal(:);
nTrials = length(coords);
nDims = size(coords{1}, 2);

% precalculate the occupancy and signal maps (unsmoothed - filter will be
% our optimization parameter),
% also linear indices of the testCoordinates
nRows = length(binEdges{1})-1;
occMaps = zeros(nRows, 1, nTrials);
fMaps = zeros(nRows, 1, nTrials);
testXYLinearIdx = cell(nTrials, 1);
for iTrial = 1:nTrials
    if isempty(coords{iTrial})
        continue;
    end
    
    [occMaps(:,:,iTrial), binCentres] = buildOccupMap(coords{iTrial}, binEdges);
    fMaps(:,:,iTrial) = buildAccumMap(coords{iTrial}, signal{iTrial}, binEdges);
    %     ind = nan(size(coords{iTrial}));
    [~, ind] = histc(coords{iTrial}, binEdges{1});
    %     [~, ind(:, 2)] = histc(coords{iTrial}(:,2), binEdges{2});
    testXYLinearIdx{iTrial} = ind;
end

occMapsXOR = zeros(size(occMaps));
fMapsXOR = zeros(size(fMaps));
for iTrial = 1:nTrials
    idx = setdiff(1:nTrials, iTrial);
    occMapsXOR(:,:,iTrial) = sum(occMaps(:,:,idx), 3);
    fMapsXOR(:,:,iTrial) = sum(fMaps(:,:,idx), 3);
end

% let's run a grid search for the optimal parameters

sz = size(occMaps);
gridValues = cell(nDims, 1);
for iDim = 1:nDims
    minValue = 1;
    maxValue = sz(iDim);
    pp = 3;
    gridValues{iDim} = linspace(minValue^(1/pp), maxValue^(1/pp), 15).^pp;
end

switch nDims
    case 1
        errValMatrix = nan(length(gridValues{1}), 1);
        
        meanVal = squeeze(sum(sum(fMapsXOR), 2)./sum(sum(occMapsXOR), 2));
        
        for mInd = 1:length(gridValues{1})
            x1 = gridValues{1}(mInd);
            errValMatrix(mInd) = mapError1D(x1, occMapsXOR, fMapsXOR, meanVal, testXYLinearIdx, signal, errPow);
        end
    otherwise
        fprintf('For now, grid search is implemented only for 2D grids\n');
end

ind = find(errValMatrix == min(errValMatrix(:)));
mOpt = ind;
stdOut = [gridValues{1}(mOpt)];
errVal = errValMatrix(ind);

% % now run the optimization
% initStd = 5*ones(1, nDims);
% [stdOut, errVal] = fminsearch(@(x) mapError(x, occMaps, fMaps, testXYLlinIdx, testSignal), initStd);

end % estOptimalStd()

%==========================================================================
function errValue = mapError1D(x1, occMaps, fMaps, meanF, testXYIdx, testF, errPow)

x1(x1<0.2) = 0.2; % HACK! This value is in pixels, not real units

hGauss = ndGaussian(x1);

[nRows, nColumns, nTrials] = size(occMaps);
chunkErrors = nan(nTrials, 1);

epsilon = 0.01;
meanFMatrix = repmat(permute(meanF, [2, 3, 1]), nRows, nColumns, 1);
% allPredictionMaps = (imfilter(fMaps, hGauss)+epsilon*meanFMatrix)./(imfilter(occMaps, hGauss)+epsilon);
fMapsSpecial = reshape(fMaps, nRows, []);
occMapsSpecial = reshape(occMaps, nRows, []);
meanFSpecial = reshape(meanFMatrix, nRows, []);
allPredictionMaps = (conv2(hGauss, 1, fMapsSpecial, 'same')+epsilon*meanFSpecial)./(conv2(hGauss, 1, occMapsSpecial, 'same')+epsilon);
allPredictionMaps = reshape(allPredictionMaps, nRows, nColumns, nTrials);

for iTrial = 1:nTrials
    if isempty(testF{iTrial})
        % there is no data for this trial
        continue;
    end
    %     occM = occMaps(:,:,iTrial);
    %     fM = fMaps(:,:,iTrial);
    %     predictionMap = filterAndDivideMaps1D(occM, fM, meanF(iTrial), hGauss);
    predictionMap = allPredictionMaps(:,:,iTrial);
    
    predictedF = predictionMap(testXYIdx{iTrial});
    
    %     chunkErrors(iChunk) = var(predictedF - testF{iChunk})/var(testF{iChunk});
    meanFTraining = meanF(iTrial);
    % the idea here is to check how much better the map is relative to just
    % using mean F map.
    switch errPow
        case 1
            chunkErrors(iTrial) = sum(abs(predictedF - testF{iTrial}))/sum(abs(meanFTraining - testF{iTrial}));
        case 2
            chunkErrors(iTrial) = sum((predictedF - testF{iTrial}).^2)/sum((meanFTraining - testF{iTrial}).^2);
        otherwise
            chunkErrors(iTrial) = sum(abs(predictedF - testF{iTrial}).^errPow)/sum(abs(meanFTraining - testF{iTrial}).^errPow);
    end
end
% errValue = nanmedian(chunkErrors);
errValue = nanmean(chunkErrors); % mean() works better, maps are smoother and nicer, "outliers have less effect"

end % mapErrorFaster()

%==========================================================================
function [stdOut, errVal, errValMatrix, gridValues] = estOptimalStd25D(coords2D, signal2D, binEdges, errPow)

% coords = coords(:);
% signal = signal(:);
[nTrials, nSubsets] = size(coords2D);
coords = cell(nTrials, 1);
signal = cell(nTrials, 1);
coordsXOR2D = cell(nTrials, nSubsets);
signalXOR2D = cell(nTrials, nSubsets);
for iTrial = 1:nTrials
    coords{iTrial, 1} = cat(1, coords2D{iTrial, 1}, coords2D{iTrial, 2});
    signal{iTrial, 1} = cat(1, signal2D{iTrial, 1}, signal2D{iTrial, 2});
    idx = setdiff(1:nTrials, iTrial);
    for iSubset = 1:nSubsets
        coordsXOR2D{iTrial, iSubset} = cell2mat(coords2D(idx, iSubset));
        signalXOR2D{iTrial, iSubset} = cell2mat(signal2D(idx, iSubset));
    end
end
nDims = size(coords{1}, 2);

% precalculate the occupancy and signal maps (unsmoothed - filter will be
% our optimization parameter),
% also linear indices of the testCoordinates
nRows = length(binEdges{1})-1;
nColumns = length(binEdges{2})-1;
occMaps = zeros(nRows, nColumns, nTrials);
fMaps = zeros(nRows, nColumns, nTrials);
testXYLinearIdx = cell(nTrials, nSubsets);
for iTrial = 1:nTrials
    if isempty(coords{iTrial})
        continue;
    end
    
    [occMaps(:,:,iTrial), binCentres] = buildOccupMap(coords{iTrial}, binEdges);
    fMaps(:,:,iTrial) = buildAccumMap(coords{iTrial}, signal{iTrial}, binEdges);
    for iSubset = 1:nSubsets
        ind = nan(size(coords2D{iTrial, iSubset}));
        [~, ind(:, 1)] = histc(coords2D{iTrial, iSubset}(:,1), binEdges{1});
        [~, ind(:, 2)] = histc(coords2D{iTrial, iSubset}(:,2), binEdges{2});
        testXYLinearIdx{iTrial, iSubset} = sub2ind([nRows, nColumns], ind(:,1), ind(:,2));
        indTrain = nan(size(coordsXOR2D{iTrial, iSubset}));
        [~, indTrain(:, 1)] = histc(coordsXOR2D{iTrial, iSubset}(:,1), binEdges{1});
        [~, indTrain(:, 2)] = histc(coordsXOR2D{iTrial, iSubset}(:,2), binEdges{2});
        trainXYLinearIdx{iTrial, iSubset} = sub2ind([nRows, nColumns], indTrain(:,1), indTrain(:,2));
    end
end
% these are training sets (XOR because they exclude trial #1 data)
occMapsXOR = zeros(size(occMaps));
fMapsXOR = zeros(size(fMaps));
for iTrial = 1:nTrials
    idx = setdiff(1:nTrials, iTrial);
    occMapsXOR(:,:,iTrial) = sum(occMaps(:,:,idx), 3);
    fMapsXOR(:,:,iTrial) = sum(fMaps(:,:,idx), 3);
    meanVal(iTrial, 1) = mean(cell2mat(signal2D(idx, 1)));
    meanVal(iTrial, 2) = mean(cell2mat(signal2D(idx, 2)));
    meanVal(iTrial, 3) = mean(cell2mat(signal(idx)));
end

% let's run a grid search for the optimal parameters

sz = size(occMaps);
gridValues = cell(nDims, 1);
for iDim = 1:nDims
    minValue = 1;
    maxValue = sz(iDim);
    pp = 3;
    gridValues{iDim} = linspace(minValue^(1/pp), maxValue^(1/pp), 15).^pp;
end

if nDims == 2
    errValMatrix = nan(length(gridValues{1}), length(gridValues{2}));
    for nInd = 1:length(gridValues{2})
        % prefilter along the 2nd dimension
        x2 = gridValues{2}(nInd);
        hGauss2 = ndGaussian(x2);
        
        occMapsXORSpecial = permute(occMapsXOR, [2, 1, 3]);
        occMapsXORSpecial = reshape(occMapsXORSpecial, nColumns, nRows*nTrials);
        occMapsXORFilt = conv2(hGauss2, 1, occMapsXORSpecial, 'same');
        occMapsXORFilt = reshape(occMapsXORFilt, nColumns, nRows, nTrials);
        occMapsXORFilt = permute(occMapsXORFilt, [2, 1, 3]);
        
        fMapsXORSpecial = permute(fMapsXOR, [2, 1, 3]);
        fMapsXORSpecial = reshape(fMapsXORSpecial, nColumns, nRows*nTrials);
        fMapsXORFilt = conv2(hGauss2, 1, fMapsXORSpecial, 'same');
        fMapsXORFilt = reshape(fMapsXORFilt, nColumns, nRows, nTrials);
        fMapsXORFilt = permute(fMapsXORFilt, [2, 1, 3]);
        
        for mInd = 1:length(gridValues{1})
            x1 = gridValues{1}(mInd);
            errValMatrix(mInd, nInd) = ...
                mapError25D(x1, occMapsXORFilt, fMapsXORFilt, meanVal, trainXYLinearIdx, signalXOR2D, testXYLinearIdx, signal2D, errPow);
        end
    end
else
    fprint('For now, grid search is implemented only for 2D grids\n');
end

ind = find(errValMatrix == min(errValMatrix(:)));
[mOpt, nOpt] = ind2sub(size(errValMatrix), ind);
stdOut = [gridValues{1}(mOpt), gridValues{2}(nOpt)];
errVal = errValMatrix(ind);

% % now run the optimization
% initStd = 5*ones(1, nDims);
% [stdOut, errVal] = fminsearch(@(x) mapError(x, occMaps, fMaps, testXYLlinIdx, testSignal), initStd);

end % estOptimalStd()

%==========================================================================
function errValue = mapError25D(x1, occMaps, fMaps, meanF, trainXYIdx, trainF, testXYIdx, testF, errPow)

x1(x1<0.2) = 0.2; % HACK! This value is in pixels, not real units

hGauss = ndGaussian(x1);

[nRows, nColumns, nTrials] = size(occMaps);
chunkErrors = nan(nTrials, 1);

epsilon = 0.01;
meanFMatrix = repmat(permute(meanF(:,3), [2, 3, 1]), nRows, nColumns, 1);
% allPredictionMaps = (imfilter(fMaps, hGauss)+epsilon*meanFMatrix)./(imfilter(occMaps, hGauss)+epsilon);
fMapsSpecial = reshape(fMaps, nRows, []);
occMapsSpecial = reshape(occMaps, nRows, []);
meanFSpecial = reshape(meanFMatrix, nRows, []);
allPredictionMaps = (conv2(hGauss, 1, fMapsSpecial, 'same')+epsilon*meanFSpecial)./(conv2(hGauss, 1, occMapsSpecial, 'same')+epsilon);
allPredictionMaps = reshape(allPredictionMaps, nRows, nColumns, nTrials);

for iTrial = 1:nTrials
    if isempty(testF{iTrial})
        % there is no data for this trial
        continue;
    end
    
    predictionMap = allPredictionMaps(:,:,iTrial);
    
    predictedF{1} = predictionMap(testXYIdx{iTrial, 1});
    predictedF{2} = predictionMap(testXYIdx{iTrial, 2});
    trainPredictedF{1} = predictionMap(trainXYIdx{iTrial, 1});
    trainPredictedF{2} = predictionMap(trainXYIdx{iTrial, 2});
    pp = getLinearScaling(trainPredictedF{1}, trainF{iTrial, 1});
    predictedF{1} = predictedF{1}*pp(1) + pp(2);
    pp = getLinearScaling(trainPredictedF{2}, trainF{iTrial, 2});
    predictedF{2} = predictedF{2}*pp(1) + pp(2);
    
    meanFTraining = meanF(iTrial);
    % the idea here is to check how much better the map is relative to just
    % using mean F map.
    switch errPow
        case 1
            chunkErrors(iTrial) = sum(abs(predictedF - testF{iTrial}))/sum(abs(meanFTraining - testF{iTrial}));
        case 2
            err1 = sum((predictedF{1}-testF{iTrial, 1}).^2);
            base1 = sum((testF{iTrial, 1}-meanF(iTrial, 1)).^2);
            err2 = sum((predictedF{2}-testF{iTrial, 2}).^2);
            base2 = sum((testF{iTrial, 2}-meanF(iTrial, 2)).^2);
            len1 = length(predictedF{1});
            len2 = length(predictedF{2});
            
            chunkErrors(iTrial) = (err1/base1*len1+err2/base2*len2)/(len1+len2);
            %             sum((predictedF - testF{iTrial}).^2)/sum((meanFTraining - testF{iTrial}).^2);
        otherwise
            chunkErrors(iTrial) = sum(abs(predictedF - testF{iTrial}).^errPow)/sum(abs(meanFTraining - testF{iTrial}).^errPow);
    end
end
% errValue = nanmedian(chunkErrors);
errValue = nanmean(chunkErrors); % mean() works better, maps are smoother and nicer, "outliers have less effect"

end % mapErrorFaster()

%============================================================================
function  out = getLinearScaling(xx, yy)

mY = mean(yy);
mX = mean(xx);
% sc = pinv(xx-mX)*(yy-mY);
sc = ((xx-mX)'*(yy-mY))/((xx-mX)'*(xx-mX));
out(1) = sc;
out(2) =  mY-sc*mX;
out = out(:);

end
