function [fMatrix] = buildZBinnedResiduals(obj, iPlane, trialIdx, zEdges)

% [~, zInd] = intersect(obj.dataTMaze.SESSION.allTrials(1).pospars, 'Z');

% pre-allocating
nTrials = length(trialIdx);
nPoints = length(zEdges)-1;
nCells = obj.nROIs(iPlane);
fMatrix = nan(nTrials, nPoints, nCells);

% nSamplesAccum = 0;
%%

data = obj.residualData{iPlane};
for trialNum = 1:nTrials
    iTrial = trialIdx(trialNum);
    if isempty(data(iTrial).residuals)
        continue;
    end
    
    zVector = data(iTrial).z;
    [occMap, ~, ~] = histcounts(zVector, zEdges);
    fRes = data(iTrial).residuals;
    for iCell=1:nCells
        binnedF = buildAccumMap(zVector, fRes(:,iCell), {zEdges})./occMap(:);
        fMatrix(trialNum, :, iCell) = binnedF';
    end
end


end 
