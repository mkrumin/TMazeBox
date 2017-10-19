function [micronsX, micronsY, micronsZ] = getPxXYZ(info)

% this function loads the header of a tiff and returns the XYZ coordinate
% of each Pixels in microns

% 2017-04-06 - LFR created

try
    allTiffInfo = dir([info.folder2pLocal, filesep, info.basename2p, '*.tif']);
    tiffName = allTiffInfo(1).name;
    filename=fullfile(info.folder2pLocal, tiffName);
    [~, header]=img.loadFrames(filename, 1, 1, 1);
catch
    fprintf('Getting the tiff from the server (local tiffs do not exist)...\n');
    allTiffInfo = dir([info.folder2p, filesep, info.basename2p, '*.tif']);
    tiffName = allTiffInfo(1).name;
    filename=fullfile(info.folder2p, tiffName);
    [~, header]=img.loadFrames(filename, 1, 1, 1);
end
% getting some parameters from the header
hh=header{1};

values = getVarFromHeader(hh, ...
    {'scanFramePeriod', 'scanZoomFactor', 'scanLinesPerFrame', 'scanPixelsPerLine'});
scanFramePeriod = str2double(values{1});
scanZoomFactor = str2double(values{2});
scanLinesPerFrame = str2double(values{3});
scanPixelsPerLine = str2double(values{4});

if isfield(info, 'microID')
    [fovx, fovy] = ppbox.zoom2fov(scanZoomFactor, infoROI.microID);
else
    [fovx, fovy] = ppbox.zoom2fov(scanZoomFactor);
end


xMicronsPerPixel = fovx / scanPixelsPerLine;
yMicronsPerPixel = fovy / scanLinesPerFrame;

micronsX = squeeze((1:scanPixelsPerLine)*xMicronsPerPixel);
micronsY = squeeze((1:scanLinesPerFrame)*yMicronsPerPixel);


if info.nPlanes>1
    micronsZ = zeros(scanLinesPerFrame, info.nPlanes);
    % only do the calculation if the piezo was moving (monty-python imaging)
    try
        % first trying to load a local copy of Timeline
        load(fullfile(info.folderTLLocal, info.basenameTL));
    catch
        load(fullfile(info.folderTL, info.basenameTL));
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
    nSamples = nCycles*info.nPlanes*scanLinesPerFrame+1;
    tInterpolated = linspace(tAxis(1), tAxis(end), nSamples);
    sInterpolated = interp1(tAxis, signal, tInterpolated);
    % averaging across cycles to get the mean cycles
    oneMeanCycle = mean(reshape(sInterpolated(1:end-1), info.nPlanes*scanLinesPerFrame, nCycles), 2);
    % extracting the bit that is corresponding the the plane of interest
    for iPlane = 1:info.nPlanes
    thisPlaneZ = oneMeanCycle((iPlane-1)*scanLinesPerFrame+1:iPlane*scanLinesPerFrame);
    
    % extracting the z coordinates of the ROI centres
        % the scaling factor is 40 um/Volt
        % the overall range is 10 Volts
     micronsZ(:, iPlane) = 40*(thisPlaneZ - referenceVoltage);
    end
else
    % nPlanes == 1, so piezo was not moving
     micronsZ = zeros(scanLinesPerFrame, info.nPlanes);
end

end




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
end