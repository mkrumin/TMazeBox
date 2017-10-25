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
