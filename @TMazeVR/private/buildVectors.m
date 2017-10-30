function [thVector, zVector, fVector, tVector] = buildVectors(obj, trialIdx, tData, fData, options)

skipT = false;
skipF = false;

if nargin<3 || isempty(tData)
    tData = [];
    tVector = [];
    skipT = true;
end
if nargin < 4 || isempty(fData)
    fData = [];
    fVector = [];
    skipF = true;
end
if nargin < 5 || ~isfield(options, 'econ')
    % if options.econ is true, then the output vectors will have the time
    % resolution of the 2p imaging, and not of the monitor refresh rate
    % false by default, for back-compatibility
    options.econ = false;
end

[~, zInd] = intersect(obj.dataTMaze.SESSION.allTrials(1).pospars, 'Z');
[~, thInd] = intersect(obj.dataTMaze.SESSION.allTrials(1).pospars, 'theta');

nSamples = 0;
for trialNum = 1:length(trialIdx)
    iTrial = trialIdx(trialNum);
    if ~isempty(obj.timesVRframes(iTrial).idx)
        if options.econ
            tFirst = obj.timesVRframes(iTrial).t(2);
            tLast = obj.timesVRframes(iTrial).t(end-1);
            nSamples = nSamples + sum(tData>=tFirst & tData <=tLast);
        else
            nSamples = nSamples + length(obj.timesVRframes(iTrial).idx)-1;
        end
    end
end

% pre-allocating
thVector = nan(nSamples, 1);
zVector = nan(nSamples, 1);
if ~skipF
    fVector = nan(nSamples, size(fData, 2));
end
if ~skipT
    tVector = nan(nSamples, 1);
end

nSamplesAccum = 0;
for trialNum = 1:length(trialIdx)
    iTrial = trialIdx(trialNum);
    tt = obj.timesVRframes(iTrial).t;
    if isempty(tt)
        continue;
    end
    idx = obj.timesVRframes(iTrial).idx(2:end);
    tt = tt(2:end-1);
    if options.econ
        tFirst = tt(1);
        tLast = tt(end);
        tIdx = tData>=tFirst & tData<=tLast;
        tAxis = tData(tIdx);
        nSamplesThisTrial = length(tAxis);
        sampleIdx = nSamplesAccum+1:nSamplesAccum+nSamplesThisTrial;
        thVector(sampleIdx) = interp1(tt, ...
            obj.dataTMaze.SESSION.allTrials(iTrial).posdata(idx, thInd), tAxis);
        zVector(sampleIdx) = interp1(tt, ...
            -obj.dataTMaze.SESSION.allTrials(iTrial).posdata(idx, zInd), tAxis);
        if ~skipF
            fVector(sampleIdx, :) = fData(tIdx, :);
        end
        if ~skipT
            tVector(sampleIdx) = tAxis(:);
        end
    else
        nSamplesThisTrial = length(tt);
        sampleIdx = nSamplesAccum+1:nSamplesAccum+nSamplesThisTrial;
        thVector(sampleIdx) = obj.dataTMaze.SESSION.allTrials(iTrial).posdata(idx, thInd);
        zVector(sampleIdx) = -obj.dataTMaze.SESSION.allTrials(iTrial).posdata(idx, zInd);
        if ~skipF
            fVector(sampleIdx, :) = interp1(tData, fData, tt);
        end
        if ~skipT
            tVector(sampleIdx) = tt;
        end
    end
    nSamplesAccum = nSamplesAccum + nSamplesThisTrial;
end

thVector = thVector*180/pi; % transform from radians to degrees

end % buildVectors()
