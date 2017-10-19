function info = registerSinglePlaneLFRdev(info, options)

% this function is doing the registration of the multiple planes of a 2p
% dataset

if nargin<2 || ~isfield(options, 'targetFrame')
    options.targetFrame='auto';
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

if isfield(options, 'cropArea')   %%LFR 27.10.14
    doCrop = true;
    area = options.cropArea(:, options.cropTag);
    area(3) = area(3) + area(1)-1;
    area(4) = area(2)+ area(4)-1;
    
    try
        filePath = fullfile(info.folderProcessed, info.basenameReg);
    catch
        info.basenameReg = sprintf('%s_plane%03d_registered', info.basename2p, options.iPlane);
        filePath = fullfile(info.folderProcessed, info.basenameReg);
        
        fprintf('registering file %s...\n', info.basenameReg);
    end
    
else
    doCrop =false;
    try
        filePath = fullfile(info.folderProcessed, info.basenameRaw);
    catch
        info.basenameRaw = sprintf('%s_plane%03d_raw', info.basename2p, options.iPlane);
        filePath = fullfile(info.folderProcessed, info.basenameRaw);
        
        fprintf('registering file %s...\n', info.basenameRaw);
    end
    
    
    
end




nFrames2Skip=[1, 0]; % number of frames to skip in the beginning/end of the stack
% the first 1-2 piezo cycles may be unstable, and the last cycle can be corrupt

[sz, ~, info] = loadArrInfo(filePath);
nFrames = sz(3);

% cutting the required frames out
indStart=nFrames2Skip(1)+1;
indEnd=nFrames-nFrames2Skip(2);
nFrames=nFrames-sum(nFrames2Skip);
% updating the planeFrame indices after removing the undesired frames
info.planeFrames=info.planeFrames(indStart:indEnd);
info.planeHeaders=info.planeHeaders(indStart:indEnd);
info.meanIntensity = info.meanIntensity(indStart:indEnd);

% this is the Gaussian filter for registration frames
hGauss = fspecial('gaussian', [5 5], 1);

% the following code is adopted (and then adapted) from regTranslations()
if isequal(options.targetFrame, 'auto')
    
    nFrames2Use = min(nFrames, 200);
    %select random frames
    frames2Use = randperm(nFrames-1, nFrames2Use)+1;
    % or evenly spaced frames
    %frames2Use = round(linspace(0.5*nFrames/nFrames2Use, (nFrames2Use-0.5)*nFrames/nFrames2Use, nFrames2Use));
    data = nan([sz(1:2), nFrames2Use]);
    for iFrame = 1:nFrames2Use
        data(:,:,iFrame) =  loadMovieFrames(filePath, frames2Use(iFrame), frames2Use(iFrame));
    end
        
    % crop when requested, LFR 27.10.14
    if doCrop        
        data = data(area(2):area(4), area(1):area(3),:);
    end
        
    
    options.targetFrame = selectTarget(1, data, nFrames2Use, 20, 1);
    
    %     %first compute a smoothed mean of each frame
    %     meanF = smooth(info.meanIntensity);
    %     %now look in the middle third of the image frames for the minimum
    %     fromFrame = round(length(info.meanIntensity)*1/3);
    %     toFrame = round(length(info.meanIntensity)*2/3);
    %     [~, idx] = min(meanF(fromFrame:toFrame));
    %     minFrame = fromFrame + idx - 1;
    %     frame = loadMovieFrames(filePath, minFrame, minFrame);
    %     %Gaussian filter the target image
    %     options.targetFrame = single(imfilter(frame, hGauss, 'same', 'replicate'));
    
elseif isnumeric(options.targetFrame) && length(options.targetFrame) == 1
    frame = loadMovieFrames(filePath, options.targetFrame, options.targetFrame);
    %Gaussian filter the target image
    options.targetFrame = single(imfilter(frame, hGauss, 'same', 'replicate'));
    
end

if options.nParThreads>1
    tmp = gcp('nocreate');
    if isempty(tmp)
        ppl=parpool(options.nParThreads);
    end
end


nChunks = ceil(nFrames/options.nFramesPerChunk);

fprintf('Calculating registration parameters...\n');
for iChunk = 1:nChunks
    
    nChars = fprintf('chunk %d/%d\n', iChunk, nChunks);
    
    frameStart = indStart + (iChunk-1)*options.nFramesPerChunk;
    frameEnd = min(frameStart + options.nFramesPerChunk - 1, indEnd);
    planeData = loadMovieFrames(filePath, frameStart, frameEnd);
    
    % crop when requested, LFR 27.10.14
    if doCrop        
        planeData = planeData(area(2):area(4), area(1):area(3),:);
    end
    
    if options.nParThreads == 1
        if options.quickReg
            [dx, dy] = img.regTranslationsMKquick(single(planeData), options.targetFrame, 'noparallel');
        else
            [dx, dy] = img.regTranslationsMK(single(planeData), options.targetFrame, 'noparallel');
        end
    else
        if options.quickReg
            [dx, dy] = img.regTranslationsMKquick(single(planeData), options.targetFrame);
        else
            [dx, dy] = img.regTranslationsMK(single(planeData), options.targetFrame);
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
for iChunk=1:nChunks
    nChars = fprintf('chunk %d/%d\n', iChunk, nChunks);
    idx = (iChunk-1)*options.nFramesPerChunk+1:...
        min(iChunk*options.nFramesPerChunk, nFrames);
    frameStart = indStart + (iChunk-1)*options.nFramesPerChunk;
    frameEnd = min(frameStart + options.nFramesPerChunk - 1, indEnd);
    planeData = loadMovieFrames(filePath, frameStart, frameEnd);
    
   if doCrop        
        planeData = planeData(area(2):area(4), area(1):area(3),:);
    end
    
    if options.translateAbs
        [mov, ~, ~]=img.translateAbs(single(planeData), info.dx(idx), info.dy(idx));
    else
        [mov, ~, ~]=img.translate(single(planeData), info.dx(idx), info.dy(idx));
    end
    if iChunk == 1
        dataPrecision = 'int16';
        
        if doCrop
            info.basenameRegisteredCrop = [info.basenameRegistered, '_cropArea_', num2str(options.cropTag)];
            fid = fopen([fullfile(info.folderProcessed,info.basenameRegisteredCrop), '.bin'], 'w');
        else
            info.basenameRegistered = strrep(info.basenameRaw, 'raw', 'registered');
            fid = fopen([fullfile(info.folderProcessed, info.basenameRegistered), '.bin'], 'w');
        end
        
        try
            fwrite(fid, int16(mov(info.validY, info.validX, :)), dataPrecision);
        catch ex
            fclose(fid);
            rethrow(ex);
        end
    else
        try
            fwrite(fid, int16(mov(info.validY, info.validX, :)), dataPrecision);
        catch ex
            fclose(fid);
            rethrow(ex);
        end
        
    end
    fprintf(repmat('\b', 1, nChars));
    
end
fclose(fid);
info.targetFrame = options.targetFrame;
s.arrSize = [length(info.validY), length(info.validX), length(info.dx)];
s.arrPrecision = dataPrecision;
if doCrop
    info.cropArea = options.cropArea;
    info.nCrop = size(options.cropArea,2);
    s.meta = info;
    save(fullfile(info.folderProcessed, info.basenameRegisteredCrop), '-struct', 's');
else
    s.meta = info;
    save(fullfile(info.folderProcessed, info.basenameRegistered), '-struct', 's');
end

fprintf('\n Finished with plane %d\n\n', options.iPlane);


% close parallel pool, if it was opened here
if options.nParThreads>1 && exist('ppl', 'var')
    delete(ppl);
end


