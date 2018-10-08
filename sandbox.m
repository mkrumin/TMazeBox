% let's construct the TMazeVR objects

[filename, folder] = uigetfile('G:\Processing\MK022', '', '', 'multiselect', 'on');
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
addpath('C:\Users\Michael\Documents\MATLAB\OldRigbox', '-begin');

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
% allExpRefs = {...
% %     '2017-09-23_1539_JL008';...
% %     '2017-09-24_1558_JL008';...
% %     '2017-09-25_1445_JL008';...
% %     '2017-09-26_1222_JL008';...
% %     '2017-09-28_1057_JL008';...
%     };

% allExpRefs = {...
%     '2017-06-12_1420_JL005';...
%     };

% allExpRefs = {...
%     '2014-08-02_2203_MK014';...
%     '2014-08-05_1937_MK014';...
%     };

allExpRefs = {...
%     '2015-07-03_1632_MK022';...
%     '2015-08-02_2051_MK022';...
%     '2015-06-19_2155_MK023';...
    '2015-07-18_154_MK023';...
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

rmpath('C:\Users\Michael\Documents\MATLAB\OldRigbox');


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
%%
%% training multiple models
% for older datasets
% addpath('\\zserver\Code\Rigging\main', '-begin');

warning off
folder = 'G:\DATA\';

% ExpRef = '2015-07-03_1632_MK022';
% ExpRef = '2015-08-02_2051_MK022';
% ExpRef = '2015-06-19_2155_MK023';
ExpRef = '2015-07-18_154_MK023';

files = dir(sprintf('G:\\DATA\\%s_TM.mat', ExpRef));
% files = dir('G:\DATA\*_TM.mat');
nFiles = length(files);
for iFile = 1:nFiles
    fprintf('iFile %2.0f/%2.0f\n', iFile, nFiles);
    load(fullfile(folder, files(iFile).name));
    ExpRef = TM.expRef;
    options.econ = true;
    for iPlane = TM.Planes
        for iROI = 1:TM.nROIs(iPlane)
            fprintf('%s Plane #%d/%d, cell #%d/%d\n', ExpRef, iPlane, length(TM.Planes), iROI, TM.nROIs(iPlane))
            TM.trainSFZDModels(iPlane, iROI, options);
        end
    end
%     fprintf('Calculating residuals..');
%     TM.getResiduals;
%     fprintf('.done\n');
    save(fullfile(folder, [ExpRef, '_TMwFits.mat']), 'TM');
end
warning on
% rmpath('\\zserver\Code\Rigging\main');

%% characterizing multiple models
warning off
folder = 'G:\DATA\';

% ExpRef = '2015-07-03_1632_MK022';
% ExpRef = '2015-08-02_2051_MK022';
% ExpRef = '2015-06-19_2155_MK023';
ExpRef = '2015-07-18_154_MK023';
files = dir(sprintf('G:\\DATA\\%s_TMwFits.mat', ExpRef));
% files = dir('G:\DATA\*_TMwFits.mat');
nFiles = length(files);
for iFile = 1:nFiles
    tic
    fprintf('iFile %2.0f/%2.0f\n', iFile, nFiles);
    load(fullfile(folder, files(iFile).name));
    ExpRef = TM.expRef;
    TM.characterizeSFZDModels;
    try
        save(fullfile(folder, [ExpRef, '_TMwExtras.mat']), 'TM');
    catch
        fprintf('Saving in ''-v7'' failed, trying to save with ''-v7.3'' flag..');
        save(fullfile(folder, [ExpRef, '_TMwExtras.mat']), 'TM', '-v7.3');
        fprintf('.done\n')
    end
    toc
end
warning on
