function [ scanPixelsPerLine , scanLinesPerFrame ] = getNumberOfPixels(info)

% this function loads the header of a tiff and returns the size of the 
% FOV in pixels

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