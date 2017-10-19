% this script is showing how to use functions included in +ppbox package to
% do basic pre-processing of the 2p data
% it was written having single-channel multiple-plane B-Scope data in mind
% adaptation to 2-channel tiffs is in the plans

% run this script section-by-section to figure out how the functions should
% be used

% by Michael Krumin, started early 2014

close all;
clear;
clc;

subject='M140422_MK008';
expDate = '2014-05-30';
exp=2;

% initializing 

info=ppbox.infoPopulate(subject, expDate, exp)


% return
%% load the datasets and divide into sub-datasets with different planes

gcp;
tic
clear options;
%options.noTiff=true;
%options.fastSave=true;
options.nFramesPerChunk = 1024;
%options.quickReg = true;
% options.translateAbs = false;
for iPlane = 3;
    options.iPlane = iPlane;
        infoRaw = ppbox.extractSinglePlane(info, options);
%     [~, ~, infoRaw] = loadArrInfo(fullfile(info.folderProcessed, sprintf('%s_plane%03d_raw', info.basename2p, iPlane)));
    infoReg = ppbox.registerSinglePlane(infoRaw, options);
end
toc
delete(gcp);

return;

%% now generate tiffs and check how good the registration was
% then, delete the raw data, probably also plane 1 and continue
% use ppbox.bin2tiff to generate tiffs

%% now register other experiments using the same targetFrames

% getting the target frames
subject='M140422_MK008';
expDate = '2014-05-27';
exp=1;
Planes = 3; % 1:5

info=ppbox.infoPopulate(subject, expDate, exp);

for iPlane = Planes
    [~, ~, infoReg{iPlane}] = loadArrInfo(fullfile(info.folderProcessed, sprintf('%s_plane%03d_registered', info.basename2p, iPlane)));
end


subject='M140422_MK008';
expDate = '2014-05-27';
exp=2;

info=ppbox.infoPopulate(subject, expDate, exp);

gcp
tic
clear options;
options.nFramesPerChunk = 1024;
options.quickReg = true;
for iPlane = Planes
    options.iPlane = iPlane;
    infoRaw = ppbox.extractSinglePlane(info, options);
    options.targetFrame = infoReg{iPlane}.targetFrame;
    infoRegistered = ppbox.registerSinglePlane(infoRaw, options);
end
toc
delete(gcp);

return;

%% now generate tiffs and check how good the registration was
% then, delete the raw data, probably also plane 1 and continue
% use ppbox.bin2tiff to generate tiffs

%% ROI detection + calcium traces extraction

subject='M140422_MK008';
expDate = '2014-05-27';
exp=2;

info=ppbox.infoPopulate(subject, expDate, exp);
CellRadius = 4; % 4 for zoom 2 experiments, 6 for zoom 3
DontPause = 0;
frames2load = 2000; % this number should be tweaked on every system (as many frames as possible, but while staying within the memory limits
for iPlane = 3:3
    basenameRegistered = sprintf('%s_plane%03d_registered', info.basename2p, iPlane);
    filePath = fullfile(info.folderProcessed, basenameRegistered);
    [sz, datatype, infoReg] = loadArrInfo(filePath);
    nFrames = sz(3);
    
%     [data, infoReg] = loadMovieFrames(filePath, 1, min(nFrames, frames2load));
%     [ROI.CellMaps, ROI.CellClasses] = AutoROIDev(data, [], CellRadius, [], [], DontPause);

    % use filename instead of loading data
    filePathMovie = fullfile(info.folderProcessed, [infoReg.chData(1).basename '_registered']);
    [ROI.CellMaps, ROI.CellClasses] = AutoROIDevMKDoNotTouch(filePathMovie, [], CellRadius, [], [], DontPause, 0, posMain, posTrace, 0.6);

    clear data;
    infoROI = infoReg;
    infoROI.ROI = ROI;
    % the following rectangles define the top-left and bottom-right corners
    % of the FOVs
    rectSource = [infoReg.validY(1), infoReg.validX(1), infoReg.validY(end), infoReg.validX(end)];
    rectDest = [1 1, size(infoReg.targetFrame)]; % full size of the targetFrame
    infoROI.targetFrameROI = remapROI(ROI.CellMaps, rectSource, rectDest);
    
    % now apply ROI chunk-by-chunk
    framesPerChunk = 2048;
    nChunks = ceil(nFrames/framesPerChunk);
    
    for iChunk = 1:nChunks
        frameFrom = framesPerChunk*(iChunk-1)+1;
        frameTo = min(nFrames, framesPerChunk*iChunk);
        [data, ~] = loadMovieFrames(filePath, frameFrom, frameTo);
        if iChunk == 1;
            infoROI.F = ppbox.applyROI(data, ROI.CellMaps);
        else
            infoROI.F = [infoROI.F; ppbox.applyROI(data, ROI.CellMaps)];
        end
        clear data;
    end
    s.arrPrecision = datatype;
    s.arrSize = sz;
    s.meta = infoROI;
    fileSavePath = strrep(filePath, 'registered', 'ROI');
%     save(fileSavePath, '-struct', 's');
    
end

%% apply the same ROIs on a different dataset - these datasets MUST have the same targetFrame

subject='M140422_MK008';
expDate = '2014-05-27';
exp=1;

oldInfo = ppbox.infoPopulate(subject, expDate, exp);

subject='M140422_MK008';
expDate = '2014-05-27';
exp=2;

% subject='MK008';
% expDate='2014-05-22';
% exp=1604;

newInfo = ppbox.infoPopulate(subject, expDate, exp);

for iPlane = 3:3
    % calculate ROI masks for other data sets (with the same targetFrame)
    oldBasename = sprintf('%s_plane%03d_ROI', oldInfo.basename2p, iPlane);
    oldFilePath = fullfile(oldInfo.folderProcessed, oldBasename);
    [~, ~, oldPlaneInfo] = loadArrInfo(oldFilePath);
    
    basename = sprintf('%s_plane%03d_registered', newInfo.basename2p, iPlane);
    filePath = fullfile(newInfo.folderProcessed, basename);
    [sz, datatype, infoReg] = loadArrInfo(filePath);

    infoROI = infoReg;
    infoROI.ROI.CellClasses = oldPlaneInfo.ROI.CellClasses;
    % the following rectangles define the top-left and bottom-right corners
    % of the FOVs
    rectSource = [oldPlaneInfo.validY(1), oldPlaneInfo.validX(1), oldPlaneInfo.validY(end), oldPlaneInfo.validX(end)]; 
    rectDest = [infoReg.validY(1), infoReg.validX(1), infoReg.validY(end), infoReg.validX(end)];
    infoROI.ROI.CellMaps = ppbox.remapROI(oldPlaneInfo.ROI.CellMaps, rectSource, rectDest);
    infoROI.targetFrameROI = oldPlaneInfo.targetFrameROI;
    
    % and then apply these ROI masks chunk-by-chunk
    nFrames = sz(3);
    framesPerChunk = 2048;
    nChunks = ceil(nFrames/framesPerChunk);
    
    for iChunk = 1:nChunks
        frameFrom = framesPerChunk*(iChunk-1)+1;
        frameTo = min(nFrames, framesPerChunk*iChunk);
        [data, ~] = loadMovieFrames(filePath, frameFrom, frameTo);
        if iChunk == 1;
            infoROI.F = ppbox.applyROI(data, infoROI.ROI.CellMaps);
        else
            infoROI.F = [infoROI.F; ppbox.applyROI(data, infoROI.ROI.CellMaps)];
        end
        clear data;
    end
    s.arrPrecision = datatype;
    s.arrSize = sz;
    s.meta = infoROI;
    fileSavePath = strrep(filePath, 'registered', 'ROI');
    save(fileSavePath, '-struct', 's');
    
end

