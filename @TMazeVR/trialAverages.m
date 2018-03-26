function [fMatrix, fModel] = trialAverages(obj, iPlane)

options.econ = true;
nTrials = obj.dataTMaze.nTrials;

% for iPlane=obj.Planes
    tData = obj.times2p{iPlane};
    fData = obj.data2p{iPlane}.F;
    f0 = prctile(fData, 20);
    fData = bsxfun(@rdivide, bsxfun(@minus, fData, f0), f0);
    nROIs = obj.nROIs(iPlane);
    fMatrix = nan(nTrials, nROIs);
    fModel = nan(nTrials, nROIs);
    for iTrial = 1:nTrials
        [~, ~, fVector, ~] = buildVectors(obj, iTrial, tData, fData, options);
        if isempty(fVector)
            continue;
        end
        fModel(iTrial, :) = mean(fVector - obj.residualData{iPlane}(iTrial).residuals, 1);
        fMatrix(iTrial, :) = mean(fVector, 1);
    end
    
%     ev{iPlane} = nan(nROIs, 1);
%     for iROI = 1:nROIs
%         tmp = corrcoef(fModel(:,iROI), fMatrix(:, iROI), 'rows', 'complete');
%         ev{iPlane}(iROI) = tmp(2).^1;
%     end
% end