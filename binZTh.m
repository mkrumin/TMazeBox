function [zAxis, thMatrix] = binZTh(z, th, nBins)

if nargin<3
    nBins = 30;
end

nTrials = length(z);

zVector = cell2mat(z);
zAxis = linspace(min(zVector), max(zVector), nBins);
dZ = diff(zAxis(1:2));
zEdges = (0:nBins)*dZ + zAxis(1) - dZ/2;
thMatrix = nan(nBins, nTrials);

for iTrial = 1:nTrials
    [N, ~, bin] = histcounts(z{iTrial}, zEdges);
    thMatrix(:, iTrial) = ...
        full(sparse(bin, ones(size(bin)), th{iTrial}, nBins, 1)./N(:));
end

