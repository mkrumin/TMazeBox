function infoUpd = redoRegistration(info, options)

% This function wil redo the registration from scratch using dx, dy
% parameters present in the info structure

% Its intended use is to be able to exactly repeat the registration done
% previously (so that we don't need to keep registered arrays after
% extracting traces)

% info - info structure with all the relevant details (see infoPopulate())
% options - an options structure

% Michael Krumin, July 2014

% 2014-07 - MK Created

%% setting the default options

if nargin<2 || ~isfield(options, 'noTiff')
    % do not save the tiff file (not really needed for unregistered data)
    options.noTiff=true;
end
if nargin<2 || ~isfield(options, 'channels')
    % by default exctract all channels
    options.channels = 1:info.nChannels;
end

%% start the processing
% get all the relevant tif-file names in that folder

% the absolute paths might have changed, we update them
infoUpd = ppbox.infoPopulate(info.subject, info.expDate, info.exp);

if isfield(infoUpd, 'folder2pLocal')
    localFiles = true;
    allTiffInfo = dir([infoUpd.folder2pLocal, filesep, info.basename2p, '*.tif']);
    if isempty(allTiffInfo)
        localFiles = false;
        fprintf('There were no locally saved tiff files at %s\n', infoUpd.folder2pLocal);
        fprintf('Will now load the files from the server (slower)\n');
        allTiffInfo = dir([infoUpd.folder2p, filesep, info.basename2p, '*.tif']);
    end
else
    localFiles = false;
    allTiffInfo = dir([infoUpd.folder2p, filesep, info.basename2p, '*.tif']);
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
    nFrames=cell(nParts, 1);
    
    iPlane = info.iPlane;
    
    if nAcqs>1
        acqString = sprintf('_acq%03d', iAcq);
    else
        acqString = '';
    end
    
    % Creating the '\Processed' sub folder, if it doesn't exist
    if ~exist(infoUpd.folderProcessed, 'dir')
        mkdir(infoUpd.folderProcessed);
    end
    
    fprintf('Extracting and reRegistering plane %d\n', iPlane);
    
    fids = cell(1, length(options.channels));
    for iCh = 1:length(options.channels)
        chID = options.channels(iCh);
        if info.nChannels>1
            chString = sprintf('_channel%03d', chID);
        else
            chString = '';
        end
        ind = strfind(info.basenameRegistered, '_rect');
        frameSetFilename = [info.basenameRegistered(1:ind-1) ...
            sprintf('_channel%03d', chID) info.basenameRegistered(ind:end)];
%         frameSetFilename = sprintf('%s%s_plane%03d%s_registered', ...
%             info.basename2p, acqString, iPlane, chString);
        fids{iCh} = fopen([fullfile(infoUpd.folderProcessed, frameSetFilename), '.bin'], 'w');
    end
    
    nFramesAccum=0;
    for iPart=1:nParts
        if localFiles
            filename=fullfile(infoUpd.folder2pLocal, tiffNames{fileIdxSorted{iAcq}(iPart)});
        else
            filename=fullfile(infoUpd.folder2p, tiffNames{fileIdxSorted{iAcq}(iPart)});
        end
        fprintf('Loading part %d/%d\n', iPart, nParts);
        if isempty(nFrames{iPart})
            if localFiles
                nFrames{iPart}=img.nFrames(fullfile(infoUpd.folder2pLocal, tiffNames{fileIdxSorted{iAcq}(iPart)}));
            else
                nFrames{iPart}=img.nFrames(fullfile(infoUpd.folder2p, tiffNames{fileIdxSorted{iAcq}(iPart)}));
            end
        end
        
        for iCh = 1:length(options.channels)
            chID = options.channels(iCh);
            
            % first frame in the current tiff file which belongs to the
            % current plane and current channel
            framesIdx = info.chData(chID).tiffFrames>nFramesAccum & ...
                info.chData(chID).tiffFrames<=nFramesAccum+nFrames{iPart};
            frames2load =  info.chData(chID).tiffFrames(framesIdx) - nFramesAccum;
            firstFrame=frames2load(1);
            lastFrame=frames2load(end);
            stride = mean(diff(frames2load));
            [data, ~]=img.loadFrames(filename, firstFrame, lastFrame, stride);
            dataPrecision = class(data);

            fprintf('Translating part %d/%d\n', iPart, nParts);
            [regData, ~, ~]=img.translate(single(data), info.dx(framesIdx), info.dy(framesIdx));
        
            try
                fprintf('Saving part %d/%d\n', iPart, nParts);
                fwrite(fids{iCh}, int16(regData(info.validY, info.validX, :)), dataPrecision);
            catch ex
                fclose(fids{iCh});
                rethrow(ex);
            end
        end
        nFramesAccum=nFramesAccum+nFrames{iPart};
    end
    for iCh = 1:length(options.channels)
        fclose(fids{iCh});
    end
    % I'm not sure we need these
    %     s.arrSize = [h w length(info.meanIntensity)];
    %     s.arrPrecision = dataPrecision;
    %     s.meta = info;
    %     save(fullfile(info.folderProcessed, planeFilename), '-struct', 's');
    fprintf('Plane %d complete.\n\n', iPlane);
end
