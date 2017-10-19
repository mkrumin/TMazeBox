function CaTraces = applyROI(frames, ROIMaps);

% caTraces = applyROI(FRAMES, ROIMaps) applies ROI masks received
% from AutoROI on the movie FRAMES and outputs an nFrames x nROIs array of
% calcium traces
% Frames - Y x X x T array
% These are raw calcium traces, NOT dF/F

% Michael Krumin
% 2014/03/19 - MK Created


nCells=length(ROIMaps);
nFrames=size(frames, 3);
frames = permute(reshape(frames, [], nFrames), [2 1]);
% nPixels = 
CaTraces=nan(nFrames, nCells);
for iCell=1:nCells
    CaTraces(:, iCell) = mean(frames(:, ROIMaps{iCell}), 2); 
end