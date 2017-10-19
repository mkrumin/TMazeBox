function names=getPlanesFromRaw(info, options)

% This function is dividing the (BScope) dataset to sub-datasets. 
% Each new sub-dataest contains a single plane
% names - names of the files generated
% info - info structure with all the relevant details (see infoPopulate())
% options - an options structure

% Michael Krumin, February 2014

% 2014-02 - MK Created

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

if nargin<2 || ~isfield(options, 'planes2extract')
    % by default extract all planes
    options.planes2extract=[1:info.nPlanes];
end


%% start the processing
% get all the relevant tif-file names in that folder

if isfield(info, 'folder2pLocal')
    allTiffInfo = dir([info.folder2pLocal, filesep, info.basename2p, '*.tif']);
    if isempty(allTiffInfo)
        fprintf('There were no locally saved tiff files at %s\n', info.folder2pLocal);
        fprintf('Will now load the files from the server (much slower)\n');
        allTiffInfo = dir([info.folder2p, filesep, info.basename2p, '*.tif']);
    end
else
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
    for iPlane=1:nPlanes
        if ~ismember(iPlane, options.planes2extract)
            continue;
        end
        fprintf('Extracting plane %d/%d\n', iPlane, nPlanes);
        nFramesAccum=0;
        for iPart=1:nParts
            filename=fullfile(info.folder2p, tiffNames{fileIdxSorted{iAcq}(iPart)});
            fprintf('Loading part %d/%d\n', iPart, nParts);
            if isempty(nFrames{iPart})
                nFrames{iPart}=img.nFrames(fullfile(info.folder2p, tiffNames{fileIdxSorted{iAcq}(iPart)}));
            end
            
            % first frame in the current tiff file which belongs to the
            % current plane
            firstFrame=mod(iPlane-mod(nFramesAccum, nPlanes), nPlanes);
            if firstFrame==0
                firstFrame=nPlanes;
            end
            frames2load=(firstFrame:nPlanes:nFrames{iPart})';
            lastFrame=frames2load(end);
            [data, headers]=img.loadFrames(filename, firstFrame, lastFrame, nPlanes);
            [h, w, nf] = size(data);
            if iPart==1
                planeData=data;
                planeFrames=frames2load;
                planeHeaders=headers;
                meanIntensity=[mean(reshape(data, h*w, nf))]';
            else
                planeData=cat(3, planeData, data);
                planeFrames=cat(1, planeFrames, frames2load+nFramesAccum);
                planeHeaders=[planeHeaders, headers];
                meanIntensity=cat(1, meanIntensity, [mean(reshape(data, h*w, nf))]');
            end
            nFramesAccum=nFramesAccum+nFrames{iPart};
        end
        info.iPlane = iPlane;
        info.planeFrames = planeFrames;
        info.planeHeaders = planeHeaders;
        info.meanIntensity = meanIntensity;

        if nAcqs>1
            names{iPlane} = sprintf('%s_acq%03d_plane%03d_raw', info.basename2p, iAcq, iPlane);
        else
            names{iPlane} = sprintf('%s_plane%03d_raw', info.basename2p, iPlane);
        end
        
        % Creating the '\Processed' sub folder, if it doesn't exist
        if ~exist(info.folderProcessed, 'dir')
            mkdir(info.folderProcessed);
        end
        
        disp('Saving the data file...')
        if options.fastSave
            % save without Matlab compression using burgbox (CB) functions
            saveArr(fullfile(info.folderProcessed, names{iPlane}), planeData, info);
        else
            % save using the standard Matlab save() function. Is
            % considerably slower (especially for >2GB files, which require '-v7.3' flag),
            % because tries to compress the data (the compression saves about 
            % 20-25% of space for typical B-Scope data)

            lastwarn('');
            save(fullfile(info.folderProcessed, names{iPlane}), 'planeData', 'info');
            if ~isempty(lastwarn)
                disp('Saving the data file to the new format...')
                save(fullfile(info.folderProcessed, names{iPlane}), 'planeData', 'info', '-v7.3');
            end
        end
        
        if ~options.noTiff
            disp('Saving the tiff...')
            saveastiff(planeData, [fullfile(info.folderProcessed, names{iPlane}), '.tif'])
        end
        fprintf('Plane %d complete.\n\n', iPlane);
    end
end
