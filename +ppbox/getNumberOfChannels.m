function [nChannels, channelColors] = getNumberOfChannels(info)

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
str = hh(strfind(hh, 'channelsSave = '):end);
ind = strfind(str, 'scanimage');
ch = str2num(str(16 : ind(1)-1));
nChannels = length(ch);

str = hh(strfind(hh, 'channelsMergeColor = '):end);
ind = strfind(str, '}');
colors = textscan(str(23:ind(1)-1), '%s', 'delimiter', ';');
if length(colors{1}) == 1
    colors = textscan(str(23:ind(1)-1), '%s', 'delimiter', ' ');
end
channelColors = cell(1, length(ch));
for k = 1:length(ch)
    channelColors{k} = colors{1}{ch(k)}(2:end-1);
end
