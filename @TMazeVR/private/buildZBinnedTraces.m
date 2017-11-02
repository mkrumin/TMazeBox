function [fMatrix] = buildZBinnedTraces(obj, trialIdx, tData, fData, zEdges)

[~, zInd] = intersect(obj.dataTMaze.SESSION.allTrials(1).pospars, 'Z');

% pre-allocating
nTrials = length(trialIdx);
nPoints = length(zEdges)-1;
nCells = size(fData, 2);
fMatrix = nan(nTrials, nPoints, nCells);

% nSamplesAccum = 0;
%%

for trialNum = 1:nTrials
    iTrial = trialIdx(trialNum);
    tt = obj.timesVRframes(iTrial).t;
    if isempty(tt)
        continue;
    end
    idx = [obj.timesVRframes(iTrial).idx(2:end); obj.timesVRframes(iTrial).idx(end)+1];
    % correction to exclude the long static frame at the end of the trial
    % (when the reward is consumed and the next trial is being prepared)
    idx = idx(1:end-1);
    tt = tt(2:end-1);
    
    zVector = -obj.dataTMaze.SESSION.allTrials(iTrial).posdata(idx, zInd);
    [occMap, ~, ~] = histcounts(zVector, zEdges);
    fInterp = interp1(tData, fData, tt);
    for iCell=1:nCells
        binnedF = buildAccumMap(zVector, fInterp(:,iCell), {zEdges})./occMap(:);
        fMatrix(trialNum, :, iCell) = binnedF';
    end
end


end % buildVectors()
