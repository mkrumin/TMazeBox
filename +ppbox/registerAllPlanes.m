function names=registerAllPlanes(info, options)

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
if nargin<2 || ~isfield(options, 'planes2register')
    % register all planes by default
    options.planes2register=[1:info.nPlanes];
end


% we will use here the fast loading from binary files
if options.fastSave
    dirData=dir(fullfile(info.folderProcessed, '*raw.bin'));
else
    dirData=dir(fullfile(info.folderProcessed, '*raw.mat'));
end
nFiles=length(dirData);
if nFiles~=info.nPlanes
    disp('Number of files is different from number of planes, check this is OK');
end
nFrames2Skip=[1, 0]; % number of frames to skip in the beginning/end of the stack
% the first 1-2 piezo cycles may be unstable, and the last cycle can be corrupt

if options.nParThreads>1
    ppl=parpool(options.nParThreads);
end

for iFile=1:nFiles

    % a hack to skip unwanted planes
    % assuming that iFile==iPlane (might not be true)
    if ~ismember(iFile, options.planes2register)
        continue;
    end
    
    filename=dirData(iFile).name;
    fprintf('loading file %s...\n', filename);
    [~, basename, ~]=fileparts(filename);
    %     load(sprintf('targetFrame_%03d.mat', iFile));
    if options.fastSave
        [planeData, meta]=loadArr(fullfile(info.folderProcessed, basename));
        % getting the plane specific info
        info=meta;
    else
        load(fullfile(info.folderProcessed, filename));
        % this command will overwrite 'planeData' and 'info'
    end
    
    [h, w, nFrames]=size(planeData);
    
    % cutting the required frames out
    indStart=nFrames2Skip(1)+1;
    indEnd=nFrames-nFrames2Skip(2);
    planeData=planeData(:,:, indStart:indEnd);
    nFrames=nFrames-sum(nFrames2Skip);
    % updating the planeFrame indices after removing the undesired frames
    info.planeFrames=info.planeFrames(indStart:indEnd);
    
    % now we will calculate the required displacements
    fprintf('Calculating registration parameters...\n');
    if options.nParThreads == 1
        [info.dx, info.dy, info.targetFrame]=...
            img.regTranslationsMK(single(planeData), options.targetFrame,...
            'noparallel', 'FramesPerChunk', num2str(options.nFramesPerChunk));
    else
        [info.dx, info.dy, info.targetFrame]=...
            img.regTranslationsMK(single(planeData), options.targetFrame, 'FramesPerChunk', num2str(options.nFramesPerChunk));
    end
    
    % now do the translations required
    fprintf('Applying registration...\n');
    nChunks=ceil(nFrames/options.nFramesPerChunk);
    for iChunk=1:nChunks
        nChars = fprintf('chunk %d/%d\n', iChunk, nChunks);
        idx=(iChunk-1)*options.nFramesPerChunk+1:min(iChunk*options.nFramesPerChunk, nFrames);
        [mov, vX, vY]=img.translate(single(planeData(:,:,idx)), info.dx(idx), info.dy(idx));
        planeData(:,:,idx)=int16(mov);
        fprintf(repmat('\b', 1, nChars));
    end
    
    % If requested, clip the frames to the maximum fully valid region
    if options.doClipping
        fprintf('Clipping the movie...\n');
        
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
    
    planeData = planeData(info.validY, info.validX, :);
    
    % save the result and return the filename
    fprintf('Saving the registered data...\n');
    
    if options.fastSave
        names{iFile}=strrep(basename, 'raw', 'registered');
        saveArr(fullfile(info.folderProcessed, names{iFile}), planeData, info);
    else
        names{iFile}=strrep(basename, 'raw', 'registered');
        lastwarn('');
        save(fullfile(info.folderProcessed, names{iFile}), 'planeData', 'info');
        if ~isempty(lastwarn)
            disp('Saving the data file to the new format...')
            save(fullfile(info.folderProcessed, names{iFile}), 'planeData', 'info', '-v7.3');
        end
    end
    
    if ~options.noTiff
        fprintf('Saving the registered tiff file...\n');
        tiffname=sprintf('%s.tiff', names{iFile});
        saveastiff(planeData, fullfile(info.folderProcessed, tiffname))
    end
    
    fprintf('Finished with plane %d\n\n', iFile);
end

% close parallel pool, if it was opened
if options.nParThreads>1
    delete(ppl);
end


