function CaTraces = applyROILFR(frames, ROIMaps)

% caTraces = applyROI(FRAMES, ROIMaps) applies ROI masks received
% from AutoROI on the movie FRAMES and outputs an nFrames x nROIs array of
% calcium traces
% Frames - Y x X x T array
% These are raw calcium traces, NOT dF/F

% Michael Krumin
% 2014/03/19 - MK Created

%LFR modified to get out pixelwise F

nCells = numel(ROIMaps);
nPx = cellfun(@numel,ROIMaps);
nFrames=size(frames, 3);
frames = permute(reshape(frames, [], nFrames), [2 1]);
% nPixels = 
CaTraces=nan(nFrames, sum(nPx));

for iCell=1:nCells
    if iCell == 1
        startPx = 1;
        endPx = nPx(iCell);
    else
        startPx = 1+nPx(iCell-1);
        endPx = nPx(iCell-1) + nPx(iCell);
    end
    CaTraces(:, startPx:endPx) = frames(:, ROIMaps{iCell}); 
end