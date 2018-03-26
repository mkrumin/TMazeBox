function compareSimilarTraces(o)

nTrials = o.dataTMaze.nTrials;
iPlane = 4;

tData = o.times2p{iPlane};
fData = o.data2p{iPlane}.F;
f0 = prctile(fData, 20);
fData = bsxfun(@minus, fData, f0);
fData = bsxfun(@rdivide, fData, f0);

options.econ = true;
resVector = cell(0);
for iTrial = 1:nTrials
    [thVector{iTrial}, zVector{iTrial}, fVector{iTrial}, ~] = ...
        o.buildVectors(iTrial, tData, fData, options);
    resVector{iTrial} = o.residualData{iPlane}(iTrial).residuals;
end
validTrials = find(~cellfun(@isempty, zVector));
idxR = o.dataTMaze.report(validTrials) == 'R';
idxL = o.dataTMaze.report(validTrials) == 'L';
thVector = thVector(validTrials);
zVector = zVector(validTrials);
fVector = fVector(validTrials);
resVector = resVector(validTrials);
nTrials = length(zVector);

zMin = max(cellfun(@min, zVector));
zMax = min(cellfun(@max, zVector));
nBins = 10;
zBinEdges = linspace(zMin - eps, zMax + eps, nBins+1);

nROIs = size(fData, 2);
zBinned = nan(nTrials, nBins);
thBinned = nan(nTrials, nBins);
fBinned = nan(nTrials, nBins, nROIs);
resBinned = nan(nTrials, nBins, nROIs);

for iTrial = 1:nTrials
    [N, ~, bin] = histcounts(zVector{iTrial}, zBinEdges);
    % TODO make sure sparse vectors are of full length
    th = sparse(1, bin+1, thVector{iTrial}, 1, nBins+1);
    thBinned(iTrial, :) = full(th(2:end))./N;
    for iROI = 1:nROIs
    f = sparse(1, bin+1, fVector{iTrial}(:, iROI), 1, nBins+1);
    fBinned(iTrial, :, iROI) = full(f(2:end))./N;
    r = sparse(1, bin+1, resVector{iTrial}(:, iROI), 1, nBins+1);
    resBinned(iTrial, :, iROI) = full(r(2:end))./N;
    end
end

%% processing zBins separately
% idxRL = bsxfun(@minus, idxR(:), idxR(:)');
% for iBin = 1:nBins
%     th = thBinned(:, iBin);
%     thDiff = bsxfun(@minus, th, th');
%     fDiff = bsxfun(@minus, fBinned(:, iBin, :), permute(fBinned(:, iBin, :), [2 1 3]));
% %     fDiff = reshape(fDiff, [], size(fDiff, 3));
%     % adding a diagonal of NaNs
%     thDiff = thDiff + 1-(1-eye(nTrials))./(1-eye(nTrials));
%     thDiff = abs(thDiff);
%     thDiffBinEdges = 0:5:60;
%     [N, ~, bin] = histcounts(thDiff, thDiffBinEdges);
%     RLDiff = nan(nTrials, nTrials, nROIs);
%     RLmeanDiff = nan(numel(thDiffBinEdges)-1, nROIs);
%     for iThDiff=1:max(bin(:))
%         ind = (idxRL.*double(bin==iThDiff));
%         RLDiff = bsxfun(@times, fDiff, ind);
%         RLDiff(repmat(~ind, 1, 1, nROIs)) = nan;
%         RLmeanDiff(iThDiff, :) = nanmean(reshape(RLDiff, [], nROIs));
%     end
%     dThAxis = (thDiffBinEdges(2:end) + thDiffBinEdges(1:end-1))/2;
%     figure
%     plot(dThAxis, RLmeanDiff(:, 4));
%     xlabel('|\Delta\theta| [deg]');
%     ylabel('<F_R - F_L>');
%     title(sprintf('z within [%2.0f %2.0f] cm', zBinEdges(iBin), zBinEdges(iBin+1)));
% end

%% let's do all the zBins together
idxRL = repmat(int8(bsxfun(@minus, idxR(:), idxR(:)')), 1, 1, nBins);
th = single(reshape(thBinned, nTrials, 1, nBins));
thDiff = bsxfun(@minus, th, permute(th, [2 1 3]));
fBinnedMat = single(reshape(fBinned, nTrials, 1, nBins, nROIs));
resBinnedMat = single(reshape(resBinned, nTrials, 1, nBins, nROIs));
fDiff = bsxfun(@minus, fBinnedMat, permute(fBinnedMat, [2 1 3 4]));
% adding a diagonal of NaNs - we don't want to compare data from 'whithin' trials
thDiff = bsxfun(@plus, thDiff, diag(nan(nTrials, 1)));
thDiff = abs(thDiff);
thDiffBinEdges = 0:5:60;
[N, ~, bin] = histcounts(thDiff, thDiffBinEdges);
nThDiffBins = numel(thDiffBinEdges)-1;
RLDiffmean = nan(nThDiffBins, nROIs);
RLDiffstd = nan(nThDiffBins, nROIs);
RLnSamples = nan(nThDiffBins, 1);
for iThDiff=1:max(bin(:))
    ind = (idxRL.*int8(bin==iThDiff));
    RLDiff = bsxfun(@times, fDiff, single(ind));
    valid = isnumeric(ind) & (ind ~= 0);
    RLDiff(repmat(~valid, 1, 1, 1, nROIs)) = nan;
    RLDiffmean(iThDiff, :) = nanmean(reshape(RLDiff, [], nROIs));
    RLDiffstd(iThDiff, :) = nanstd(reshape(RLDiff, [], nROIs));
    RLnSamples(iThDiff) = sum(valid(:))/2;
end
%% now repeat for the residuals
resRLDiffmean = nan(nThDiffBins, nROIs);
resRLDiffstd = nan(nThDiffBins, nROIs);
resRLnSamples = nan(nThDiffBins, 1);
resDiff = bsxfun(@minus, resBinnedMat, permute(resBinnedMat, [2 1 3 4]));
for iThDiff=1:max(bin(:))
    ind = (idxRL.*int8(bin==iThDiff));
    RLDiff = bsxfun(@times, resDiff, single(ind));
    valid = isnumeric(ind) & (ind ~= 0);
    RLDiff(repmat(~valid, 1, 1, 1, nROIs)) = nan;
    resRLDiffmean(iThDiff, :) = nanmean(reshape(RLDiff, [], nROIs));
    resRLDiffstd(iThDiff, :) = nanstd(reshape(RLDiff, [], nROIs));
    resRLnSamples(iThDiff) = sum(valid(:))/2;
end
    
%%
dThAxis = (thDiffBinEdges(2:end) + thDiffBinEdges(1:end-1))/2;

EV = 1-[o.trainingData{iPlane}.errVals]';
[~, idx] = sort(EV, 'descend');

for iCell = 1:10
    iROI = idx(iCell);
    figure
    subplot(1, 3, 1);
    zAxis = o.trainingData{iPlane}(iROI).zThetaBinCentres{1};
    thAxis = o.trainingData{iPlane}(iROI).zThetaBinCentres{2};
    imagesc(thAxis, zAxis, o.trainingData{iPlane}(iROI).zThetaMap)
    axis xy equal tight 
    title('The z-\theta Model');

    ax = subplot(1, 3, 2);
    %     plot(dThAxis, RLDiffmean(:, 4));
    errorbar(dThAxis, RLDiffmean(:, iROI), RLDiffstd(:, iROI)./sqrt(RLnSamples));
    hold on;
    errorbar(dThAxis, resRLDiffmean(:, iROI), resRLDiffstd(:, iROI)./sqrt(resRLnSamples));
    ax.XLim = ([min(thDiffBinEdges), max(thDiffBinEdges)]);
    ax.YLim = ([-max(abs(ax.YLim)), max(abs(ax.YLim))]);
    plot(ax.XLim, [0 0], 'k:')
    xlabel('|\Delta\theta| [deg]');
    ylabel('<F_R - F_L>');
    title(sprintf('iROI = %2.0f', iROI));
    legend('Data', 'Residuals');
    
    subplot(1, 3, 3);
    imagesc(thAxis, zAxis, o.trainingData{iPlane}(iROI).residualMap)
    caxis([min( o.trainingData{iPlane}(iROI).zThetaMap(:)), ...
        max( o.trainingData{iPlane}(iROI).zThetaMap(:))])
    axis xy equal tight
    title('z-\theta map of residuals');

end
%%
keyboard;

%%
figure
errorbar(zBinEdges(2:end), mean(fBinned(idxR, :, iROI))', std(fBinned(idxR, :, iROI))');
hold on;
errorbar(zBinEdges(2:end), mean(fBinned(idxL, :, iROI))', std(fBinned(idxL, :, iROI))');