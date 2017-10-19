function info = extractSinglePlaneZ(info, options)

% This function is dividing the (BScope) dataset to sub-datasets.
% Each new sub-dataset contains a single plane
% info - info structure with all the relevant details (see infoPopulate())
% options - an options structure

% Michael Krumin, February 2014

% 2014-02 - MK Created (getPlanesFromRaw)
% 2014-06 - MK modified for low mem (or large arrays) use (extractSinglePlane)

%% setting the default options

if nargin<2 || ~isfield(options, 'noTiff')
    % do not save the tiff file (not really needed for unregistered data)
    options.noTiff=true;
end

if nargin<2 || ~isfield(options, 'fastSave')
    % set to true if you want to use the fast array saving instead of
    % matlab save() function. Will be much faster, but the data is not
    % compressed (you loose approx 20-25% of space for typical BScope data)
    % make sure you have saveArr and loadArr (writeen by Chris) if you use
    % this fastSave option
    options.fastSave=true;
end

if nargin<2 || ~isfield(options, 'iPlane')
    % by default extract the first plane
    options.iPlane=1;
end

if nargin<2 || ~isfield(options, 'channels')
    % by default exctract all channels
    options.channels = 1:info.nChannels;
end


%% start the processing
% get all the relevant tif-file names in that folder

if isfield(info, 'folder2pLocal')
    localFiles = true;
    allTiffInfo = dir([info.folder2pLocal, filesep, info.basename2p, '*.tif']);
    if isempty(allTiffInfo)
        localFiles = false;
        fprintf('There were no locally saved tiff files at %s\n', info.folder2pLocal);
        fprintf('Will now load the files from the server (slower)\n');
        allTiffInfo = dir([info.folder2p, filesep, info.basename2p, '*.tif']);
    end
else
    localFiles = false;
    allTiffInfo = dir([info.folder2p, filesep, info.basename2p, '*.tif']);
end

nFiles = length(allTiffInfo);
tiffNames = cell(nFiles, 1);
acqNumber = nan(nFiles, 1);
partNumber = nan(nFiles, 1);
for iFile = 1:nFiles
    tiffNames{iFile} = allTiffInfo(iFile).name;
    acqNumber(iFile) = str2num(tiffNames{iFile}(end-10:end-8));
    partNumber(iFile) = str2num(tiffNames{iFile}(end-6:end-4));
end

acquisitions=unique(acqNumber);
nAcqs=length(acquisitions);

if nAcqs>1
    warning('There is more than one acquisition in this folder, something might be wrong');
end

for iAcq=1:nAcqs
    [partsSorted{iAcq} fileIdxSorted{iAcq}]=sort(partNumber(acqNumber==acquisitions(iAcq)));
end

for iAcq=1:nAcqs
    
    nParts=length(fileIdxSorted{iAcq});
    nPlanes=info.nPlanes;
    nFrames=cell(nParts, 1);
    nFrameSets = (info.nPlanes-1) * info.nChannels;
    
    iPlane = options.iPlane;
    info.iPlane = iPlane;
    
    
    if nAcqs>1
        acqString = sprintf('_acq%03d', iAcq);
    else
        acqString = '';
    end
    planeFilename = sprintf('%s%s_plane%03d', ...
        info.basename2p, acqString, iPlane);
    
    % Creating the '\Processed' sub folder, if it doesn't exist
    if ~exist(info.folderProcessed, 'dir')
        mkdir(info.folderProcessed);
    end
    
    fprintf('Extracting plane %d/%d\n', iPlane, nPlanes);
    
    fids = cell(1, length(options.channels));
    for iCh = 1:length(options.channels)
        chID = options.channels(iCh);
        if info.nChannels>1
            chString = sprintf('_channel%03d', chID);
        else
            chString = '';
        end
        frameSetFilename = sprintf('%s%s_plane%03d%s_raw', ...
            info.basename2p, acqString, iPlane, chString);
        basename = sprintf('%s%s_plane%03d%s', ...
            info.basename2p, acqString, iPlane, chString);
        fids{iCh} = fopen([fullfile(info.folderProcessed, frameSetFilename), '.bin'], 'w');
        info.chData(chID).basename = basename;
    end
    
    nFramesAccum=0;
    for iPart=1:nParts
        if localFiles
            filename=fullfile(info.folder2pLocal, tiffNames{fileIdxSorted{iAcq}(iPart)});
        else
            filename=fullfile(info.folder2p, tiffNames{fileIdxSorted{iAcq}(iPart)});
        end
        fprintf('Loading part %d/%d\n', iPart, nParts);
        if isempty(nFrames{iPart})
            if localFiles
                nFrames{iPart}=img.nFrames(fullfile(info.folder2pLocal, tiffNames{fileIdxSorted{iAcq}(iPart)}));
            else
                nFrames{iPart}=img.nFrames(fullfile(info.folder2p, tiffNames{fileIdxSorted{iAcq}(iPart)}));
            end
        end
        
        for iCh = 1:length(options.channels)
            chID = options.channels(iCh);
    
            % first frame in the current tiff file which belongs to the
            % current plane and current channel
            iFrameSet = (iPlane-1) * info.nChannels + chID;
                        
            firstFrame= mod(iFrameSet-mod(nFramesAccum, nFrameSets), nFrameSets);
        
            if firstFrame==0
                firstFrame=nFrameSets;
            end
                        
            frames2load=(firstFrame:nFrameSets:nFrames{iPart})';
            
%             % HACK 16.02.2015
%             frames2load(frames2load > 1021) = [];
%             %
            
            lastFrame=frames2load(end);
            
            [data, headers]=img.loadFrames(filename, firstFrame, lastFrame, nFrameSets);
            [h, w, nf] = size(data);
            if iPart==1
                if iCh == 1
                    info.planeHeaders = headers(1);
                end
                dataPrecision = class(data);
                info.chData(chID).tiffFrames = frames2load;
                info.chData(chID).meanIntensity = mean(reshape(data, h*w, nf))';

                try
                    fwrite(fids{iCh}, data, dataPrecision);
                catch ex
                    fclose(fids{iCh});
                    rethrow(ex);
                end
            else
                info.chData(chID).tiffFrames = cat(1, info.chData(chID).tiffFrames, frames2load+nFramesAccum);
                info.chData(chID).meanIntensity = cat(1, info.chData(chID).meanIntensity, mean(reshape(data, h*w, nf))');
                
                try
                    fwrite(fids{iCh}, data, dataPrecision);
                catch ex
                    fclose(fids{iCh});
                    rethrow(ex);
                end
            end
        end
        nFramesAccum=nFramesAccum+nFrames{iPart};
    end
    for iCh = 1:length(options.channels)
        fclose(fids{iCh});
    end
    % for backward compatibility
    greenCh = strcmp({info.chData.color}, 'green');
    if sum(greenCh)>0 %****LFR added if clause 2.7.15
        info.meanIntensity = info.chData(greenCh).meanIntensity;
    else
        redCh = strcmp({info.chData.color}, 'red');
        info.meanIntensity = info.chData(redCh).meanIntensity;       
        
    end
    info.basenameRaw = [planeFilename '_raw'];
    
    info.basenamePlane = planeFilename;
    info.planeFrames = ceil(info.chData(options.channels(1)).tiffFrames / info.nChannels);
    s.arrSize = [h w length(info.chData(options.channels(1)).meanIntensity)];
    s.arrPrecision = dataPrecision;
    s.meta = info;
    save(fullfile(info.folderProcessed, [planeFilename '_raw']), '-struct', 's');
    fprintf('Plane %d complete.\n\n', iPlane);
end
