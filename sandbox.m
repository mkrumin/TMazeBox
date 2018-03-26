% let's construct the TMazeVR objects

[filename, folder] = uigetfile('G:\Processing\JL008\2017-07-27\', '', '', 'multiselect', 'on');
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
% for older datasets
% addpath('\\zserver\Code\Rigging\main', '-begin');

warning off
folder = 'G:\DATA\';

files = dir('G:\DATA\TM_*');
nFiles = length(files);
for iFile = 1:nFiles
    fprintf('iFile %2.0f/%2.0f\n', iFile, nFiles);
    load(fullfile(folder, files(iFile).name));
    TM = TM.TMVR;
    ExpRef = TM.expRef;
    options.econ = true;
    for iPlane = TM.Planes
        for iROI = 1:TM.nROIs(iPlane)
%             fprintf('%s Plane #%d/%d, cell #%d/%d\n', ExpRef, iPlane, length(TM.Planes), iROI, TM.nROIs(iPlane))
            TM.trainMap_CV(iPlane, iROI, options);
        end
    end
    fprintf('Calculating residuals..');
    TM.getResiduals;
    fprintf('.done\n');
    save(fullfile(folder, [ExpRef, '_TM.mat']), 'TM');
end
warning on
% rmpath('\\zserver\Code\Rigging\main');



%% training
% for older datasets
addpath('\\zserver\Code\Rigging\main', '-begin');

folder = 'G:\DATA\';

% allExpRefs = {'2017-07-15_1708_JL008';...
%     '2017-07-27_1433_JL008';...
%     '2017-07-28_1359_JL008';...
%     '2017-08-12_1056_JL008';...
%     '2017-08-14_1414_JL008';...
%     '2017-09-18_1707_JL008';...
%     '2017-09-19_1117_JL008';...
%     };

% allExpRefs = {};

% someting is wrong with these datasets, photodiode signal?
allExpRefs = {...
%     '2017-09-23_1539_JL008';...
%     '2017-09-24_1558_JL008';...
%     '2017-09-25_1445_JL008';...
%     '2017-09-26_1222_JL008';...
    '2017-09-28_1057_JL008';...
    };


for iExpRef = 1:length(allExpRefs)
    ExpRef = allExpRefs{iExpRef}
    TM = TMazeVR(ExpRef);
    options.econ = true;
    for iPlane = TM.Planes
        for iROI = 1:TM.nROIs(iPlane)
            fprintf('%s Plane #%d/%d, cell #%d/%d', ExpRef, iPlane, length(TM.Planes), iROI, TM.nROIs(iPlane))
            TM.trainMap_CV(iPlane, iROI, options);
            fprintf(': ExplVar = %6.4f\n', 1-TM.trainingData{iPlane}(iROI).errVals);
        end
    end
    fprintf('Calculating residuals..');
    TM.getResiduals;
    fprintf('.done\n');
    save(fullfile(folder, [ExpRef, '_TM.mat']), 'TM');
end

rmpath('\\zserver\Code\Rigging\main');


%% loading some data

allExpRefs = {'2017-07-15_1708_JL008';...
    '2017-07-27_1433_JL008';...
    '2017-07-28_1359_JL008';...
    '2017-08-12_1056_JL008';...
    '2017-08-14_1414_JL008'};

folder = 'G:\DATA\';
ExpRef = allExpRefs{5};
load(fullfile(folder, [ExpRef, '_TM.mat']))
TM.rocAnalysis;

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

%% resildual analysis
% if ~exist('TM', 'var')
%     if strfind(hostname, 'zenbook')
%     [filename, folder] = uigetfile('C:\Processing\JL008\2017-07-15\1708\*_TM.mat', 'Select TM file', '');
%     elseif strfind(hostname, 'zero')
%     [filename, folder] = uigetfile('G:\Processing\JL008\2017-07-15\1708\*_TM.mat', 'Select TM file', '');
%     end
%     load(fullfile(folder, filename));
% end

% TM.getResiduals;

%% plotting maps with the corresponding residuals

warning('off', 'MATLAB:nargchk:deprecated');
iCell = 0;
nRows = 4;
nColumns = 6;
nCellsPerFigure = nRows*nColumns/2;
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
sameCAxis = true;
for iCell = 1:3*nCellsPerFigure
    if (mod(iCell, nCellsPerFigure)==1)
        figure;
    end
    iPlane = planeIdx(ind(iCell));
    iROI = roiIdx(ind(iCell));
    thAxis = TM.trainingData{iPlane}(iROI).zThetaBinCentres{2};
    zAxis = TM.trainingData{iPlane}(iROI).zThetaBinCentres{1};
    EV = 1-errVal(ind(iCell));

    
    iPlot = 2*(mod(iCell-1, nCellsPerFigure)+1)-1;
    subplot(nRows, nColumns, iPlot);
    imagesc(thAxis, zAxis, TM.trainingData{iPlane}(iROI).rawDataMap);
    axis equal tight xy
    if sameCAxis
%         clim = minmax([TM.trainingData{iPlane}(iROI).zThetaMap(:)', ...
%             TM.trainingData{iPlane}(iROI).residualMap(:)']);
        clim = minmax([TM.trainingData{iPlane}(iROI).rawDataMap(:)', ...
            TM.trainingData{iPlane}(iROI).rawResMap(:)']);
        caxis(clim);
    else
        colorbar;
    end
    title(sprintf('Cell %1.0f, EV %4.2f', iCell, EV));
    
    iPlot = 2*(mod(iCell-1, nCellsPerFigure)+1);
    subplot(nRows, nColumns, iPlot);
    imagesc(thAxis, zAxis, TM.trainingData{iPlane}(iROI).rawResMap);
    axis equal tight xy
    if sameCAxis
        caxis(clim);
    else
        colorbar;
    end

end


