function times=getFrameTimes(info, frames)
info;
if nargin<2
    try
        frames=info.planeFrames;
    catch
        frames='all'
    end
end

try
    % first trying to load a local copy of Timeline
    load(fullfile(info.folderTLLocal, info.basenameTL));
catch
    load(fullfile(info.folderTL, info.basenameTL));
end

nInputs=length(Timeline.hw.inputs);
for iInput=1:nInputs
    if isequal(Timeline.hw.inputs(iInput).name, 'neuralFrames')
        ind=iInput;
        break;
    end
end

nTotalFrames=Timeline.rawDAQData(end, ind);
TTLs=[0; diff(Timeline.rawDAQData(:, ind))];
idx=find(TTLs);

if isequal(frames, 'all')
    frames=1:nTotalFrames;
else
    frames=frames(frames<=nTotalFrames);
end

times=Timeline.rawDAQTimestamps(idx(frames));

