function infoOut = getROIExtras(infoROI)

% this function will add the following information about the ROIs:
% (x, y, z) location of the cell
% infoOut.ROI.CellXYPixels will contain [x, y] coordinates in pixels [1, 1]
% being the top left corner of the frame
% infoOut.ROI.CellXYZMicrons will contain [x, y, z] coordinates in microns.
% x and y are relative to the top-left corner of the frame and z go
% more negative as we go deeper into the brain. Z is NOT the actual depth
% inside the brain, it is just a relative coordinate inside the 400um range
% of the piezo. MK - Seems to be buggy.
% MK Fix 2015/03/09 - now z is relative to the STARTING position of the piezo.
% This is the position where piezo sits before starting the volumetric
% scanning, and this is the position where we usually take the reference
% image (in focusing mode). This info should be enough, together with measuring the depth of
% this starting position relative to the brain surface, to infer the actual
% depth of the ROI.

% tShift - sampling time of the cell, relative to the frame onset
% infoOut.ROI.CellTShift - time in seconds that take the scanning laser beam to get
% from top of the frame to the relevant y-position. This gives an estimate
% of the more precise sampling timestamps of the ROIs.

% 2014-07 Michael Krumin

% 2016-02 F. Rossi added the call to zoom2fov to get correct estimation of
% FOV size at different zoomfactors

infoOut = infoROI;
%% get x,y coords, both in pixels and in microns

if isfield(infoROI.chData, 'targetFrame')
    [h, w] = size(infoROI.chData(infoROI.registrationChannel).targetFrame);
else
    [h, w] = size(infoROI.targetFrame);
end

values = getVarFromHeader(infoROI.planeHeaders{1}, ...
    {'scanFramePeriod', 'scanZoomFactor', 'scanLinesPerFrame', 'scanPixelsPerLine'});
scanFramePeriod = str2double(values{1});
scanZoomFactor = str2double(values{2});
scanLinesPerFrame = str2double(values{3});
scanPixelsPerLine = str2double(values{4});

if isfield(infoROI, 'microID')
    [fovx, fovy] = ppbox.zoom2fov(scanZoomFactor, infoROI.microID);
else
    [fovx, fovy] = ppbox.zoom2fov(scanZoomFactor);
end

for iROI = 1:length(infoROI.targetFrameROI)
    [ySub, xSub] = ind2sub([h, w], infoROI.targetFrameROI{iROI});
    xCentrePixels = mean(xSub);
    yCentrePixels = mean(ySub);
    %   for zoomFactor 2 we have 500um field of view, scale accordingly
    %     xMicronsPerPixel = 1/scanPixelsPerLine/scanZoomFactor*1000;
    %     yMicronsPerPixel = 1/scanLinesPerFrame/scanZoomFactor*1000;
    xMicronsPerPixel = fovx / scanPixelsPerLine;
    yMicronsPerPixel = fovy / scanLinesPerFrame;
    xCentreMicrons = xCentrePixels*xMicronsPerPixel;
    yCentreMicrons = yCentrePixels*yMicronsPerPixel;
    infoOut.ROI.CellXYPixels{iROI} = [xCentrePixels, yCentrePixels];
    infoOut.ROI.CellXYZMicrons{iROI} = [xCentreMicrons, yCentreMicrons, NaN];
end

%% get z coords ( only microns)
if infoROI.nPlanes>1
    % only do the calculation if the piezo was moving (monty-python imaging)
    try
        % first trying to load a local copy of Timeline
        load(fullfile(infoROI.folderTLLocal, infoROI.basenameTL));
    catch
        load(fullfile(infoROI.folderTL, infoROI.basenameTL));
    end
    
    nInputs=length(Timeline.hw.inputs);
    for iInput=1:nInputs
        if isequal(Timeline.hw.inputs(iInput).name, 'piezoPosition')
            indPosition=iInput;
        end
        if isequal(Timeline.hw.inputs(iInput).name, 'piezoCommand')
            indCommand=iInput;
        end
        if isequal(Timeline.hw.inputs(iInput).name, 'neuralFrames')
            indFrames=iInput;
        end
    end
    
    % here comes a bit of cumbersome code, but it seems that it works ok
    deltas = [0; -diff(Timeline.rawDAQData(:, indCommand))];
    cycleStarts = (deltas>0.7*max(deltas));
    % skip the first 5 cycles and then take ten cycles to average across
    nCycles = 10;
    nStart = 5;
    indices = find(cycleStarts, nStart + nCycles, 'first');
    if max(abs(diff(diff(indices(nStart+1:end)))))>1
        % take only reliable sequence,
        % if not reliable skip a few cycles and try again
        nStart = nStart + 5;
        indices = find(cycleStarts, nStart + nCycles, 'first');
    end
    % defining start and end of the nCycle segment
    startFrame = indices(nStart);
    endFrame = indices(end);
    % extracting the segment
    signal = Timeline.rawDAQData(startFrame:endFrame, indPosition);
    tAxis = Timeline.rawDAQTimestamps(startFrame:endFrame);
    nSamples = round(Timeline.hw.daqSampleRate/2);
    referenceVoltage = mean(Timeline.rawDAQData(1:nSamples, indPosition));
    % interpolating to fit the number of lines in the image
    nSamples = nCycles*infoROI.nPlanes*scanLinesPerFrame+1;
    tInterpolated = linspace(tAxis(1), tAxis(end), nSamples);
    sInterpolated = interp1(tAxis, signal, tInterpolated);
    % averaging across cycles to get the mean cycles
    oneMeanCycle = mean(reshape(sInterpolated(1:end-1), infoROI.nPlanes*scanLinesPerFrame, nCycles), 2);
    % extracting the bit that is corresponding the the plane of interest
    thisPlaneZ = oneMeanCycle((infoROI.iPlane-1)*scanLinesPerFrame+1:infoROI.iPlane*scanLinesPerFrame);
    
    % extracting the z coordinates of the ROI centres
    for iROI = 1:length(infoROI.targetFrameROI)
        % the scaling factor is 40 um/Volt
        % the overall range is 10 Volts
        infoOut.ROI.CellXYZMicrons{iROI}(3) = 40*(referenceVoltage - interp1(1:scanLinesPerFrame, thisPlaneZ, infoOut.ROI.CellXYPixels{iROI}(2)));
    end
    
else
    % nPlanes == 1, so piezo was not moving
    for iROI = 1:length(infoROI.targetFrameROI)
        infoOut.ROI.CellXYZMicrons{iROI}(3) = 0;
    end
end


%% get tShift

for iROI = 1:length(infoROI.targetFrameROI)
    infoOut.ROI.CellTShift{iROI} = scanFramePeriod*infoOut.ROI.CellXYPixels{iROI}(2)/scanLinesPerFrame;
end


%% this sub function will extract the values from the ScanImage tiff header
function values = getVarFromHeader(str, fields)

% str is the header
% fields is a cell array of strings with variable names
% values is a cell array of corresponding values, they will be strings

ff = strsplit(str, {' = ', 'scanimage.SI4.'});
if ~iscell(fields)
    fields = cell(fields);
end
values = cell(size(fields));

for iField = 1:length(fields)
    ind = find(ismember(ff, fields{iField}));
    values{iField} = ff{ind+1};
end