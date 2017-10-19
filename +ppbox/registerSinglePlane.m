function info = registerSinglePlane(info, options)

% this function is doing the registration of the multiple planes of a 2p
% dataset

if nargin<2 || ~isfield(options, 'targetFrame')
    options.targetFrame='auto';
end
if nargin<2 || ~isfield(options, 'nFrames4TargetSelection')
    options.nFrames4TargetSelection=100;
end
if nargin<2 || ~isfield(options, 'nFramesPerChunk')
    options.nFramesPerChunk=512;
end
if nargin<2 || ~isfield(options, 'doClipping')
    options.doClipping=true;
end
if nargin<2 || ~isfield(options, 'nParThreads')
    % define this >1 if you want to use parrallel processing
    options.nParThreads=4;
end
if nargin<2 || ~isfield(options, 'fastSave')
    % true if using fast binary saving/loading
    options.fastSave=true;
end
if nargin<2 || ~isfield(options, 'noTiff')
    % set to true if do not want to save tiffs
    options.noTiff=false;
end
if nargin<2 || ~isfield(options, 'iPlane')
    % register all planes by default
    options.iPlane=info.iPlane;
end
if nargin<2 || ~isfield(options, 'quickReg')
    options.quickReg = false;
end
if nargin<2 || ~isfield(options, 'translateAbs')
    options.translateAbs = false;
end
if nargin<2 || ~isfield(options, 'registrationChannel')
    % register red channel by default (if it exists), otherwise the first
    % analysed channel
    regCh = find(strcmp({info.chData.color}, 'red'));
    tmp = 0;
    while isempty(regCh)
        tmp = tmp + 1;
        if ~isempty(info.chData(tmp).tiffFrames)
            regCh = tmp;
        end
    end
    options.registrationChannel = regCh;
end

regCh = options.registrationChannel;
if info.nChannels>1
    chString = sprintf('_channel%03d', regCh);
else
    chString = '';
end
try
    filePath = fullfile(info.folderProcessed, [info.chData(regCh).basename '_raw']);
catch
    info.chData(regCh).basename = sprintf('%s_plane%03d%s', ...
        info.basename2p, options.iPlane, chString);
    filePath = fullfile(info.folderProcessed, [info.chData(regCh).basename '_raw']);
    
    % for backward compatibility
    greenCh = strcmp({info.chData.color}, 'green');
    if info.nChannels>1
        chString = sprintf('_channel%03d', greenCh);
    else
        chString = '';
    end
    info.basenameRaw = sprintf('%s_plane%03d%s', ...
        info.basename2p, options.iPlane, chString);
end

nFrames2Skip=[1, 0]; % number of frames to skip in the beginning/end of the stack
% the first 1-2 piezo cycles may be unstable, and the last cycle can be corrupt

fprintf('registering file %s_raw...\n', info.chData(regCh).basename);
[sz, prec, info] = loadArrInfo(fullfile(info.folderProcessed, [info.basenamePlane '_raw']));
nFrames = sz(3);

m = memmapfile([filePath, '.bin']);
m.Format =  {'int16', sz, 'frames'};

% cutting the required frames out
indStart=nFrames2Skip(1)+1;
indEnd=nFrames-nFrames2Skip(2);
nFrames=nFrames-sum(nFrames2Skip);
% updating the planeFrame indices and tiffFrame indices after removing the undesired frames
info.planeFrames=info.planeFrames(indStart:indEnd);

% for backward compatibility
info.meanIntensity = info.meanIntensity(indStart:indEnd);



for iCh = 1:info.nChannels
    if isempty(info.chData(iCh).tiffFrames)
        continue
    end
    info.chData(iCh).meanIntensity = info.chData(iCh).meanIntensity(indStart:indEnd);
    info.chData(iCh).tiffFrames = info.chData(iCh).tiffFrames(indStart:indEnd);
    
end


% this is the Gaussian filter for registration frames
hGauss = fspecial('gaussian', [5 5], 1);

% the following code is adopted (and then adapted) from regTranslations()
if isequal(options.targetFrame, 'auto')
    
    nFrames2Use = min(nFrames, options.nFrames4TargetSelection);
    % select random frames
    % frames2Use = randperm(nFrames-1, nFrames2Use)+1;
    % or evenly spaced frames
    frames2Use = round(linspace(0.5*nFrames/nFrames2Use, ...
        (nFrames2Use-0.5)*nFrames/nFrames2Use, nFrames2Use));
    
    data = m.Data.frames(:,:,frames2Use);
    clear m
    
    options.targetFrame = selectTarget(1, data, nFrames2Use, 20, 1, hGauss);
    clear data
    
elseif isnumeric(options.targetFrame) && length(options.targetFrame) == 1
    
    data = m.Data.frames(:,:,options.targetFrame);
    clear m
    %Gaussian filter the target image
    options.targetFrame = single(imfilter(data, hGauss, 'same', 'replicate'));
    
elseif isequal(options.targetFrame, 'average') %% added by LFR on 2.7.15 to deal with zstacks easily
    
    data = mean(m.Data.frames,3);
   
    clear m
    %Gaussian filter the target image
    options.targetFrame = single(imfilter(data, hGauss, 'same', 'replicate'));
    
end

if options.nParThreads>1
    try
        tmp = gcp('nocreate');
        if isempty(tmp)
            ppl=parpool(options.nParThreads);
        end
    catch
        % for older vesions of matlab (prior to 2013b)
        tmp = matlabpool('size');
        if ~tmp
            ppl = [];
            matlabpool;
        end
    end
end


nChunks = ceil(nFrames/options.nFramesPerChunk);

fprintf('Calculating registration parameters...\n');

for iChunk = 1:nChunks
    
    nChars = fprintf('chunk %d/%d\n', iChunk, nChunks);
    
    frameStart = indStart + (iChunk-1)*options.nFramesPerChunk;
    frameEnd = min(frameStart + options.nFramesPerChunk - 1, indEnd);
    planeData = loadMovieFrames(filePath, frameStart, frameEnd, sz, prec);
    
    if options.nParThreads == 1
        if options.quickReg
            [dx, dy] = img.regTranslationsMKquick(single(planeData), ...
                options.targetFrame, 'noparallel');
        else
            [dx, dy] = img.regTranslationsMK(single(planeData), ...
                options.targetFrame, 'noparallel');
        end
    else
        if options.quickReg
            [dx, dy] = img.regTranslationsMKquick(single(planeData), ...
                options.targetFrame);
        else
            [dx, dy] = img.regTranslationsMK(single(planeData), ...
                options.targetFrame);
        end
    end
    
    if iChunk == 1
        info.dx = dx(:);
        info.dy = dy(:);
    else
        info.dx = [info.dx; dx(:)];
        info.dy = [info.dy; dy(:)];
    end
    
    fprintf(repmat('\b', 1, nChars));
end

% If requested, clip the frames to the maximum fully valid region
[h, w, ~] = size(planeData);
if options.doClipping
    dxMax = max(0, ceil(max(info.dx)));
    dxMin = min(0, floor(min(info.dx)));
    dyMax = max(0, ceil(max(info.dy)));
    dyMin = min(0, floor(min(info.dy)));
    info.validX = (1 + dxMax):(w + dxMin);
    info.validY = (1 + dyMax):(h + dyMin);
else
    info.validX = 1:w;
    info.validY = 1:h;
end

% now do the translations required
fprintf('Applying registration and saving...\n');

fids = cell(1, info.nChannels);
for iCh = 1:info.nChannels
    if isempty(info.chData(iCh).tiffFrames)
        continue
    end
    fids{iCh} = fopen(fullfile(info.folderProcessed, [info.chData(iCh).basename '_registered.bin']), 'w');
end

dataPrecision = 'int16';
for iCh = 1:info.nChannels
    if isempty(info.chData(iCh).tiffFrames)
        continue
    end
    filePath = fullfile(info.folderProcessed, [info.chData(iCh).basename '_raw']);
    for iChunk=1:nChunks
        nChars = fprintf('chunk %d/%d\n', iChunk, nChunks);
        idx = (iChunk-1)*options.nFramesPerChunk+1:...
            min(iChunk*options.nFramesPerChunk, nFrames);
        frameStart = indStart + (iChunk-1)*options.nFramesPerChunk;
        frameEnd = min(frameStart + options.nFramesPerChunk - 1, indEnd);
        planeData = loadMovieFrames(filePath, frameStart, frameEnd, sz, prec);
        if options.translateAbs
            [mov, ~, ~]=img.translateAbs(single(planeData), info.dx(idx), info.dy(idx));
        else
            [mov, ~, ~]=img.translate(single(planeData), info.dx(idx), info.dy(idx));
        end
        try
            fwrite(fids{iCh}, int16(mov(info.validY, info.validX, :)), dataPrecision);
        catch ex
            fclose(fids{iCh});
            rethrow(ex);
        end
        fprintf(repmat('\b', 1, nChars));
    end
    fclose(fids{iCh});
end
info.chData(regCh).targetFrame = options.targetFrame;
info.registrationChannel = regCh;

% for backward compatibility
info.basenameRegistered = [info.basenamePlane 'registered'];
info.targetFrame = options.targetFrame;

s.arrSize = [length(info.validY), length(info.validX), length(info.dx)];
s.arrPrecision = dataPrecision;
s.meta = info;
save(fullfile(info.folderProcessed, [info.basenamePlane '_registered']), '-struct', 's');

fprintf('Finished with plane %d\n\n', options.iPlane);

% close parallel pool, if it was opened here
if options.nParThreads>1 && exist('ppl', 'var')
    if isempty(ppl)
        % for older versions of matlab
        matlabpool close;
    else
        delete(ppl);
    end
    
end


