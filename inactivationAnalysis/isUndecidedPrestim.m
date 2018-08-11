function [isUndecided, thCorrected] = isUndecidedPrestim(data, zLaserStart)

nTrials = length(data.z);
thLaserStart = nan(nTrials, 1);

for iTrial = 1:nTrials

    idx = find(data.z{iTrial}>zLaserStart, 1, 'first') - 1;
    if ~isempty(idx)
        thLaserStart(iTrial) =  data.theta{iTrial}(idx);
    end
   
end


bounds = prctile(thLaserStart, [25,75]);
med = nanmedian(thLaserStart);
thCorrected = thLaserStart - med;
% 
% mTh = nanmean(th10)
% stdTh = nanstd(th10);
% bounds(1) = mTh - (0.5*stdTh);
% bounds(2) = mTh + (0.5*stdTh);

isUndecided = (thLaserStart>bounds(1) & thLaserStart<bounds(2));