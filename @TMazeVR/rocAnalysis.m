function rocAnalysis(obj)

%% Getting the trial indices

% keyboard
trialsAll = 1:length(obj.dataTMaze.contrastSequence);

cc = unique(obj.dataTMaze.contrastSequence);
nContrasts = length(cc);
for iCC = 1:nContrasts
    trialsGoR{iCC} = find(obj.dataTMaze.contrastSequence' == cc(iCC) & obj.dataTMaze.report == 'R');
    trialsGoL{iCC} = find(obj.dataTMaze.contrastSequence' == cc(iCC) & obj.dataTMaze.report == 'L');
    trialsC{iCC} = find(obj.dataTMaze.contrastSequence' == cc(iCC) & obj.dataTMaze.outcome == 'C');
    trialsW{iCC} = find(obj.dataTMaze.contrastSequence' == cc(iCC) & obj.dataTMaze.outcome == 'W');
    trialsStimR{iCC} = find(obj.dataTMaze.contrastSequence' == cc(iCC) & obj.dataTMaze.contrastSequence'>0);
    trialsStimL{iCC} = find(obj.dataTMaze.contrastSequence' == cc(iCC) & obj.dataTMaze.contrastSequence'<0);

end
% trialsGoR{length(cc)+1} = find(obj.dataTMaze.report == 'R');
% trialsGoL{length(cc)+1} = find(obj.dataTMaze.report == 'L');

[trialsGoR, groupLabels] = groupTrials(trialsGoR, cc);
trialsGoL = groupTrials(trialsGoL, cc);
trialsC = groupTrials(trialsC, cc);
trialsW = groupTrials(trialsW, cc);
trialsStimR = groupTrials(trialsStimR, cc);
trialsStimL = groupTrials(trialsStimL, cc);

nGroups = length(trialsGoR);

%% Cutting F traces (trial-wise) for all the ROIs
fprintf('Analyzing experiment %s\n', obj.expRef);
zEdges = obj.trainingData{min(obj.Planes)}(1).zEdges;
zEdges = linspace(zEdges(1), zEdges(end), 11);
nZ = length(zEdges)-1;
fMatrix = cell(max(obj.Planes), 1);
for iPlane = obj.Planes
    tData = obj.times2p{iPlane}';
    nCells = obj.nROIs(iPlane);
    fMatrix{iPlane} = nan(length(trialsAll), nZ, nCells);
    fData = obj.data2p{iPlane}.F;
    f0 = prctile(fData, 20);
    fData = bsxfun(@rdivide, bsxfun(@minus, fData, f0), f0);
    [fMatrix{iPlane}] = buildZBinnedTraces(obj, trialsAll, tData, fData, zEdges);
end

% reshape fMatrix to be nTrials x nZBins x nCells matrix
fMatrix = cell2mat(reshape(fMatrix, 1, 1, []));

% keyboard;

%% Performing the decision-style analysis

nCells = size(fMatrix, 3);

nChars = 0;
for iCell = 1:nCells
    fprintf(repmat('\b', 1, nChars));
    nChars = fprintf('ROI %d/%d\n', iCell, nCells);
    %   take traces of one ROI
    data = fMatrix(:,:,iCell);
    % normalize data to be between 0 and 1
    minF = min(data(:));
    maxF = max(data(:));
    data = (data - minF)/(maxF-minF);
    if all(all(isnan(data)))
        % this ROI traces were not made in the previous section
        continue;
    end
    
    for iGroup = 1:nGroups
        
        % calculate the trials' statistics
        tracesR = data(trialsGoR{iGroup}, :);
        tracesL = data(trialsGoL{iGroup}, :);
        
        tracesC = data(trialsC{iGroup}, :);
        tracesW = data(trialsW{iGroup}, :);
        
        tracesStimR = data(trialsStimR{iGroup}, :);
        tracesStimL = data(trialsStimL{iGroup}, :);

        
%         % U-test
%         try
%         [pValRL(iCell, iGroup), hRL(iCell, iGroup), stat] = ranksum(nanmean(tracesR, 2), nanmean(tracesL, 2));
%         [pValCW(iCell, iGroup), hCW(iCell, iGroup), stat] = ranksum(nanmean(tracesC, 2), nanmean(tracesW, 2));
%         end
        
        % ROC analysis
        for iBin = 1:nZ
            valuesR = tracesR(:, iBin);
            valuesL = tracesL(:, iBin);
            valuesR = valuesR(~isnan(valuesR));
            valuesL = valuesL(~isnan(valuesL));
            [rocRL(iBin, iCell, iGroup), tprRL{iBin, iCell, iGroup}, fprRL{iBin, iCell, iGroup}] ...
                = rocArea(valuesR, valuesL);
            
            valuesR = tracesStimR(:, iBin);
            valuesL = tracesStimL(:, iBin);
            valuesR = valuesR(~isnan(valuesR));
            valuesL = valuesL(~isnan(valuesL));
            [rocStimRL(iBin, iCell, iGroup), tprStimRL{iBin, iCell, iGroup}, fprStimRL{iBin, iCell, iGroup}] ...
                = rocArea(valuesR, valuesL);
            
            valuesC = tracesC(:, iBin);
            valuesW = tracesW(:, iBin);
            valuesC = valuesC(~isnan(valuesC));
            valuesW = valuesW(~isnan(valuesW));
            [rocCW(iBin, iCell, iGroup), tprCW{iBin, iCell, iGroup}, fprCW{iBin, iCell, iGroup}] ...
                = rocArea(valuesC, valuesW);
        end
    end
    
end
fprintf('\n');

%% plotting is done here

nRows = 1;
nColumns = 3;

subplot(1, 3, [1 2]);
rocRL(rocRL==0) = NaN;
rocRL = abs(rocRL-0.5)+0.5;
nCells = sum(~all(isnan(rocRL), 2));
% errorbar(1:nGroups, nanmean(rocRL), nanstd(rocRL, [], 1)/sqrt(nCells))
errorbar(1:nGroups, nanmean(rocRL), nanstd(rocRL, [], 1), 'b')
hold on;
rocCW(rocCW==0) = NaN;
rocCW = abs(rocCW-0.5)+0.5;
nCells = sum(~all(isnan(rocCW), 2));
% errorbar(1:nGroups, nanmean(rocCW), nanstd(rocCW, [], 1)/sqrt(nCells))
errorbar(1:nGroups, nanmean(rocCW), nanstd(rocCW, [], 1), 'r:')

rocStimRL(rocStimRL==0) = NaN;
rocStimRL = abs(rocStimRL-0.5)+0.5;
nCells = sum(~all(isnan(rocStimRL), 2));
% errorbar(1:nGroups, nanmean(rocStimRL), nanstd(rocStimRL, [], 1)/sqrt(nCells))
errorbar(1:nGroups, nanmean(rocStimRL), nanstd(rocStimRL, [], 1), 'g--')



xlim([0.5, nGroups+0.5]);
ylim([0, 1]);
plot(xlim, [0.5, 0.5], 'k:');
set(gca, 'XTick', [1:nGroups], 'XTickLabel', groupLabels)
title(sprintf('ROC   (%d cells)', nCells));
legend('ROC_{R-L} + std', 'ROC_{C-W} + std', 'ROC_{Stim R-L} + std');
xlabel('Contrast [%]');

subplot(1, 3, 3);
showPC(obj.dataTMaze, gca);
expRefStr = strrep(obj.expRef, '_', '\_');
title(sprintf('%s', expRefStr));

set(gcf, 'name', obj.expRef);
drawnow;
% pause;

%%
function [trialsOut, groupLabels] = groupTrials(trialsIn, cc)

groupLabels = {'High % L', 'Low % L', '0 %', 'Low % R', 'High % R', 'All %'};
bins = {[-50, -25];...
    [-12, -6];...
    [0];...
    [6, 12];...
    [25, 50];...
    [cc]};


% groupLabels = {'-50%', '-25%', '-12%', '-6%', '0%', '6%', '12%', '25%', '50%', 'All %'};
% bins = { -50; -25; -12; -6; 0; 6; 12; 25; 50; cc};


% groupLabels = {'0%', '\pm6%', '\pm12%', '\pm25%', '\pm50%', 'All %'};
% bins = {[0]; [-6, 6]; [-12, 12]; [-25, 25]; [-50, 50]; cc};

trialsIn = trialsIn(:)'; % making sure it is a row vector of cell (for cell2mat later)
nGroups = length(groupLabels);

trialsOut = cell(nGroups, 1);

for iGroup = 1:nGroups
    
    [~, ccInd, ~] = intersect(cc, bins{iGroup});
    trialsOut{iGroup} = unique(cell2mat(trialsIn(ccInd)));
    
end

