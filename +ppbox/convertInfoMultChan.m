function infoMultChan = convertInfoMultChan(infoSingleChan)

infoMultChan = infoSingleChan;

infoMultChan.nChannels = 1;
infoMultChan.chData.color = 'green';
infoMultChan.registrationChannel = 1;
if isfield(infoSingleChan, 'basenameRegistered')
    infoMultChan.chData.basename = strrep(infoSingleChan.basenameRegistered, ...
        '_registered', '');
    infoMultChan.basenamePlane = strrep(infoSingleChan.basenameRegistered, ...
        '_registered', '');
end
if isfield(infoSingleChan, 'planeFrames')
    infoMultChan.chData.tiffFrames = infoSingleChan.planeFrames;
end
if isfield(infoSingleChan, 'meanIntensity')
    infoMultChan.chData.meanIntensity = infoSingleChan.meanIntensity;
end
if isfield(infoSingleChan, 'targetFrame')
    infoMultChan.chData.targetFrame = infoSingleChan.targetFrame;
end
if isfield(infoSingleChan, 'planeHeaders')
    infoMultChan.planeHeaders = infoSingleChan.planeHeaders(1);
end