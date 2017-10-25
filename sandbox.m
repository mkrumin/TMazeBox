% let's construct the TMazeVR objects

[filename, folder] = uigetfile('G:\Processing\JL008\2017-07-15\1708\', '', '', 'multiselect', 'on');
if ~iscell(filename)
    filename = {filename};
end

nFiles = length(filename);

for iFile = 1:nFiles
    data = load(fullfile(folder, filename{iFile}));
end

info = data.meta;
TM = TMazeVR(info.expRef);

%% training
for iPlane = TM.Planes
    for iROI = 1:TM.nROIs(iPlane)
        fprintf('Plane #%d/%d, cell #%d/%d\n', iPlane, length(TM.Planes), iROI, TM.nROIs(iPlane))
        TM.trainMap_CV(iPlane, iROI);
    end
end

%% plotting
warning('off', 'MATLAB:nargchk:deprecated');
iCell = 0;
nRows = 4;
nColumns = 6;
nCellsPerFigure = nRows*nColumns;
planeIdx = [];
roiIdx = [];
errVal = [];
for iPlane = TM.Planes
    planeIdx = cat(1, planeIdx, iPlane*ones(TM.nROIs(iPlane), 1));
    roiIdx = cat(1, roiIdx, [1:TM.nROIs(iPlane)]');
    errVal = cat(1, errVal, [TM.trainingData{iPlane}(:).errVals]');
end

[~, ind] = sort(errVal, 'ascend');
nCells = length(errVal);
for iCell = 1:nCellsPerFigure
    if (mod(iCell, nCellsPerFigure)==1)
        figure;
    end
    subplot(nRows, nColumns, mod(iCell-1, nCellsPerFigure)+1);
    iPlane = planeIdx(ind(iCell));
    iROI = roiIdx(ind(iCell));
    thAxis = TM.trainingData{iPlane}(iROI).zThetaBinCentres{2};
    zAxis = TM.trainingData{iPlane}(iROI).zThetaBinCentres{1};
    imagesc(thAxis, zAxis, TM.trainingData{iPlane}(iROI).zThetaMap);
    axis equal tight xy
    EV = 1-errVal(ind(iCell));
    title(EV);
end

