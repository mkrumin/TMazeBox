function info = registerSinglePlaneSmartDevMK(info, options)

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

try
    filePath = fullfile(info.folderProcessed, info.basenameRaw);
catch
    info.basenameRaw = sprintf('%s_plane%03d_raw', info.basename2p, options.iPlane);
    filePath = fullfile(info.folderProcessed, info.basenameRaw);
end

nFrames2Skip=[1, 0]; % number of frames to skip in the beginning/end of the stack
% the first 1-2 piezo cycles may be unstable, and the last cycle can be corrupt

fprintf('registering file %s...\n', info.basenameRaw);
[sz, ~, info] = loadArrInfo(filePath);
nFrames = sz(3);

m = memmapfile([filePath, '.bin']);
m.Format =  {'int16', sz, 'frames'};

% cutting the required frames out
indStart=nFrames2Skip(1)+1;
indEnd=nFrames-nFrames2Skip(2);
nFrames=nFrames-sum(nFrames2Skip);
% updating the planeFrame indices after removing the undesired frames
info.planeFrames=info.planeFrames(indStart:indEnd);
% info.planeHeaders=info.planeHeaders(indStart:indEnd);
info.meanIntensity = info.meanIntensity(indStart:indEnd);

% this is the Gaussian filter for registration frames
hGauss = fspecial('gaussian', [7 7], 1);

% the following code is adopted (and then adapted) from regTranslations()
if isequal(options.targetFrame, 'auto')
    
    nFrames2Use = min(nFrames, options.nFrames4TargetSelection);
    % select random frames
    % frames2Use = randperm(nFrames-1, nFrames2Use)+1;
    % or evenly spaced frames
    frames2Use = round(linspace(0.5*nFrames/nFrames2Use, (nFrames2Use-0.5)*nFrames/nFrames2Use, nFrames2Use));
    %     data = nan([sz(1:2), nFrames2Use]);
    %     for iFrame = 1:nFrames2Use
    %         data(:,:,iFrame) =  loadMovieFrames(filePath, frames2Use(iFrame), frames2Use(iFrame));
    %     end
    
    data = m.Data.frames(:,:,frames2Use);
    clear m;
    
    sz = size(data);
    options.targetFrame = selectTargetDev(1, data, nFrames2Use, floor(nFrames2Use/20), 1, hGauss);
    %     clear data;

    figure
    imagesc(options.targetFrame);
    axis equal tight;
    colormap gray
    hold on;
    
%     nY = 3;
%     nX = 3;
%     yy = round(linspace(0, sz(1), nY+1));
%     xx = round(linspace(0, sz(2), nX+1));
%     iRect = 0;
%     for iY = 1:nY
%         for iX = 1:nX
%             iRect = iRect+1;
%             posBank(iRect).ymin = yy(iY)+1;
%             posBank(iRect).ymax = yy(iY+1);
%             posBank(iRect).xmin = xx(iX)+1;
%             posBank(iRect).xmax = xx(iX+1);
%         end
%     end
%     
    for iRect = 1:100 % 100 looks like a reasonable maximum number of 
        h = imrect(gca);
        pos = wait(h);
        set(h, 'Visible', 'off');
        hold on;
        if (prod(pos(3:4))==0)
            %             if the size of the rectangle is 0 - exit the loop
            drawnow; 
            break;
        else
            posBank(iRect).ymin = max(1, floor(pos(2)));
            posBank(iRect).ymax = min(sz(2), ceil(pos(2)+pos(4)));
            posBank(iRect).xmin = max(1, floor(pos(1)));
            posBank(iRect).xmax = min(sz(2), ceil(pos(1)+pos(3)));
            
            %     plot(cumsum([pos(1) pos(3)]), cumsum([pos(2) pos(4)]), '.r')
            plot([pos(1) pos(1) pos(1)+pos(3) pos(1)+pos(3) pos(1)], [pos(2) pos(2)+pos(4) pos(2)+pos(4) pos(2) pos(2)], 'r:');
            yInd = round(pos(2)):round(pos(2)+pos(4));
            xInd = round(pos(1)):round(pos(1)+pos(3));
            drawnow;
        end
        
    end
    
    nRects = length(posBank);
    
    for iRect = 1:nRects
        p = posBank(iRect);
        plot([p.xmin, p.xmax, p.xmax, p.xmin, p.xmin], [p.ymin, p.ymin, p.ymax, p.ymax, p.ymin], 'g:')
    end
    drawnow;
    
    collage = zeros(sz(1), sz(2));
    targetBank = cell(nRects, 1);
    for iRect = 1:nRects
        p = posBank(iRect);
        targetBank{iRect} = ...
            selectTargetDev(1, data(p.ymin:p.ymax, p.xmin:p.xmax, :), nFrames2Use, floor(nFrames2Use/20), 1, hGauss);
        collage(p.ymin:p.ymax, p.xmin:p.xmax) = targetBank{iRect};
    end
    
    figure
    imagesc(collage);
    axis equal tight;
    colormap gray;
    drawnow;
    
elseif isnumeric(options.targetFrame) && length(options.targetFrame) == 1
    frame = loadMovieFrames(filePath, options.targetFrame, options.targetFrame);
    %Gaussian filter the target image
    options.targetFrame = single(imfilter(frame, hGauss, 'same', 'replicate'));
    
end

% pause;

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
allDx = nan(nFrames, nRects);
allDy = nan(nFrames, nRects);
for iChunk = 1:nChunks
    
    nChars = fprintf('chunk %d/%d\n', iChunk, nChunks);
    idx = (iChunk-1)*options.nFramesPerChunk+1:...
        min(iChunk*options.nFramesPerChunk, nFrames);
    
    frameStart = indStart + (iChunk-1)*options.nFramesPerChunk;
    frameEnd = min(frameStart + options.nFramesPerChunk - 1, indEnd);
    planeData = loadMovieFrames(filePath, frameStart, frameEnd);
    %     planeData = m.Data.frames(:,:,frameStart:frameEnd);
    
    for iRect = 1:nRects
        yInd = posBank(iRect).ymin:posBank(iRect).ymax;
        xInd = posBank(iRect).xmin:posBank(iRect).xmax;
        planeDataCropped = planeData(yInd, xInd, :);
        options.targetFrameCropped = targetBank{iRect};
        
        if options.nParThreads == 1
            if options.quickReg
                [dx, dy] = img.regTranslationsMKquick(single(planeDataCropped), options.targetFrameCropped, 'noparallel');
            else
                [dx, dy] = img.regTranslationsMK(single(planeDataCropped), options.targetFrameCropped, 'noparallel');
            end
        else
            if options.quickReg
                [dx, dy] = img.regTranslationsMKquick(single(planeDataCropped), options.targetFrameCropped);
            else
                [dx, dy] = img.regTranslationsMK(single(planeDataCropped), options.targetFrameCropped);
            end
        end
        
        allDx(idx, iRect) = dx(:);
        allDy(idx, iRect) = dy(:);
    end
    fprintf(repmat('\b', 1, nChars));
end

% If requested, clip the frames to the maximum fully valid region
dxMax = max(0, ceil(max(allDx)));
dxMin = min(0, floor(min(allDx)));
dyMax = max(0, ceil(max(allDy)));
dyMin = min(0, floor(min(allDy)));
% [h, w, ~] = size(planeData);
if options.doClipping
    for iRect = 1:nRects
        w = posBank(iRect).xmax - posBank(iRect).xmin + 1;
        h = posBank(iRect).ymax - posBank(iRect).ymin + 1;
        validX{iRect} = (1 + dxMax(iRect)):(w + dxMin(iRect));
        validY{iRect} = (1 + dyMax(iRect)):(h + dyMin(iRect));
    end
else
    for iRect = 1:nRects
        w = posBank(iRect).xmax - posBank(iRect).xmin + 1;
        h = posBank(iRect).ymax - posBank(iRect).ymin + 1;
        validX{iRect} = 1:w;
        validY{iRect} = 1:h;
    end
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
    %     planeData = m.Data.frames(:,:,frameStart:frameEnd);
    
    for iRect = 1:nRects
        yInd = posBank(iRect).ymin:posBank(iRect).ymax;
        xInd = posBank(iRect).xmin:posBank(iRect).xmax;
        planeDataRect = planeData(yInd, xInd, :);

        if options.translateAbs
            [mov, ~, ~]=img.translateAbs(single(planeDataRect), allDx(idx, iRect), allDy(idx, iRect));
        else
            [mov, ~, ~]=img.translate(single(planeDataRect), allDx(idx, iRect), allDy(idx, iRect));
        end
        if iChunk == 1
            dataPrecision = 'int16';
            basenameRegistered = sprintf('%s_rect%d', strrep(info.basenameRaw, 'raw', 'registered'), iRect);
            fid(iRect) = fopen([fullfile(info.folderProcessed, basenameRegistered), '.bin'], 'w');
            try
                fwrite(fid(iRect), int16(mov(validY{iRect}, validX{iRect}, :)), dataPrecision);
            catch ex
                fclose(fid(iRect));
                rethrow(ex);
            end
        else
            try
                fwrite(fid(iRect), int16(mov(validY{iRect}, validX{iRect}, :)), dataPrecision);
            catch ex
                fclose(fid(iRect));
                rethrow(ex);
            end
            
        end
    end
    fprintf(repmat('\b', 1, nChars));
    
end
% closing the data files and also saving the additional information in the .mat files
for iRect = 1:nRects
    fclose(fid(iRect));
    info.targetFrame = options.targetFrame;
    info.targetFrameCrop = targetBank{iRect};
    info.cropPosition = posBank(iRect);
    info.collage = collage;
    info.validX = validX{iRect} + posBank(iRect).xmin - 1;
    info.validY = validY{iRect} + posBank(iRect).ymin - 1;
    info.dx = allDx(:, iRect);
    info.dy = allDy(:, iRect);
    info.basenameRegistered = sprintf('%s_rect%d', strrep(info.basenameRaw, 'raw', 'registered'), iRect);
    s.arrSize = [length(info.validY), length(info.validX), length(info.dx)];
    s.arrPrecision = dataPrecision;
    s.meta = info;
    save(fullfile(info.folderProcessed, info.basenameRegistered), '-struct', 's');
end


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


