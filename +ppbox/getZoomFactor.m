function zoomFactor = getZoomFactor(info)

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
str = hh(strfind(hh, 'scanZoomFactor = '):end);
ind = strfind(str, 'scanimage');
zoomFactor = str2double(str(18 : ind(1)-1));