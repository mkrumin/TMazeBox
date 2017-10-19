function bin2tiff(fnIn, sz, precision)

% This function will convert a binary *.bin file into tiff (if the data
% inside the .bin is really a tiff image)

% You will need to interactively select .bin file(s) to convert, all the
% rest is done automatically
% Alternatively, you can supply a single (full) filename fnIn as an argument

% Tiff files will be divided into parts to keep them under 2GB, for files
% above 2GB savetiff works extremely slowly, we might consider even changing
% this to 1GB

% 2014-05 Michael Krumin


% choose file(s) interactively if no filename was supplied
if nargin<1
    DialogTitle = 'Selct file(s) to convert to .tiff';
    [FileName,PathName,FilterIndex] = uigetfile('*.bin', DialogTitle, ppbox.getTmpFolder,...
        'MultiSelect', 'on');
else
    [PathName, fn, fe] = fileparts(fnIn);
    FileName = [fn, fe];
end

if isstr(FileName)
    FileName = {FileName};
end

for iFile = 1:length(FileName)
    [~, fn, ~] = fileparts(FileName{iFile});
    fprintf('converting file %s...\n', fn);
    if nargin < 2
        [sz, precision] = loadArrInfo(fullfile(PathName, fn));
    end
    nFrames = sz(3);
    info = dir(fullfile(PathName, FileName{iFile}));
    nParts = ceil(info.bytes/2^30); % let's keep files to less than 1GB
    framesPerPart = ceil(nFrames/nParts);
    
    if ~exist('C:\Temp\', 'dir')
        mkdir('C:\Temp\');
    end
    
    currentFolder = cd('C:\Temp\'); % the saveastiff works much faster if the current folder is local and not on the network
    for iPart = 1:nParts
        frameStart = (iPart-1)*framesPerPart+1;
        frameEnd = min(iPart*framesPerPart, nFrames);
%         tic
        fprintf('loading part %d/%d..', iPart, nParts);
        data = loadMovieFrames(fullfile(PathName, fn), frameStart, frameEnd, ...
            sz, precision);
        fprintf('.done\n');
%         toc
        % only add an iPart index if there is more than one part
%         tic
        fprintf('tiffing part %d/%d..', iPart, nParts);
        if nParts>1
            saveastiff(data, fullfile(PathName, [fn, '_', num2str(iPart, '%03d'),'.tiff']))
        else
            saveastiff(data, fullfile(PathName, [fn,'.tiff']))
        end
        fprintf('.done\n');
%         toc
    end
    cd(currentFolder);
    clear data;
end
