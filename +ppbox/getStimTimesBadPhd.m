function times=getStimTimesBadPhd(info, badPhd)

%% LFR variant of getStimTimes to use SyncEcho when the photodiode signal is bad

try
    % first trying to load a local copy of Timeline
    load(fullfile(info.folderTLLocal, info.basenameTL));
catch
    load(fullfile(info.folderTL, info.basenameTL));
end

nInputs=length(Timeline.hw.inputs);

if badPhd
    for iInput=1:nInputs
        if isequal(Timeline.hw.inputs(iInput).name, 'syncEcho')
            ind=iInput;
            break;
        end
    end
else
    for iInput=1:nInputs
        if isequal(Timeline.hw.inputs(iInput).name, 'photoDiode')
            ind=iInput;
            break;
        end
    end
end

phd=Timeline.rawDAQData(:, ind); 
phd=(phd-min(phd))/(max(phd)-min(phd));
thr=0.5; % using one threshold here

above=phd>thr;
deltas=[0; diff(above)];
goingUpTimes=Timeline.rawDAQTimestamps(deltas==1);
goingDownTimes=Timeline.rawDAQTimestamps(deltas==-1);

if badPhd
goingUpTimes= goingUpTimes + 0.017; % according to MK
goingDownTimes= goingDownTimes +0.017;
end

startIdx=[];
endIdx=[];
for iUDP=1:Timeline.mpepUDPCount
    if isequal(Timeline.mpepUDPEvents{iUDP}(1:9), 'StimStart')
        startIdx=[startIdx; iUDP];
    end
    if isequal(Timeline.mpepUDPEvents{iUDP}(1:7), 'StimEnd')
        endIdx=[endIdx; iUDP];
    end
end

% the first UP after the StimStart UDP

times.onset=[];
for iUDP=1:length(startIdx)
    tmp=min(goingUpTimes(goingUpTimes>Timeline.mpepUDPTimes(startIdx(iUDP))));
    times.onset=[times.onset; tmp];
end
% the last DOWN before the StimEnd UDP
times.offset=[];
delayConst = 0.05; % a delay constant, for the case when the 'StimEnd' udp arrived before the stimulus finished to play
for iUDP=1:length(endIdx)
    tmp=max(goingDownTimes(goingDownTimes<(Timeline.mpepUDPTimes(endIdx(iUDP)) + delayConst)));
    times.offset=[times.offset; tmp];
end

% figure
% plot(Timeline.rawDAQTimestamps, phd);
% hold on;
% plot(times.onset, 0.5, 'r.');
% plot(times.offset, 0.5, 'g.');
