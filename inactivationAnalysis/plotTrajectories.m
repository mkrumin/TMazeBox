function plotTrajectories(data)

figure;

nTrials = length(data.z);
th10 = nan(nTrials, 1);
for iTrial = 1:nTrials

    plot(data.theta{iTrial}, data.z{iTrial}, 'b')
    hold on;
    
    idx = find(data.z{iTrial}>10, 1, 'first') - 1;
    if ~isempty(idx)
    th10(iTrial) =  data.theta{iTrial}(idx);
    end
end
axis tight equal

figure;
histogram(th10, 20);

med = nanmedian(th10)
ave = nanmean(th10)

%% plotting the trajectories

% prc = [25 50 75];
% zLims = [5 60];
% cc = unique(abs(contrast));
% idxFR = finished & random;
% idxCorrect = outcome == 'C';
% idxR = behavior == 'R';
% idxL = behavior == 'L';
% idxStimL = optiStim(:, 1) > 0 & optiStim(:, 2) == 0;
% idxStimR = optiStim(:, 2) > 0 & optiStim(:, 1) == 0;
% idxStimBoth = optiStim(:, 1) > 0 & optiStim(:, 2) > 0;
% idxStimNone = optiStim(:, 1) == 0 & optiStim(:, 2) == 0;
% [zAxis, thMatrix] = binZTh(allZ, allTh, 30);
% 
% figure('Name', list(1).animalName)
% 
% nRows = length(cc)-1;
% for iC = 2:length(cc)
%     idxC = (abs(contrast) == cc(iC));
%     
%     subplot(nRows, 4, 1+(iC-2)*4)
%     idxBlue = find(idxC & idxCorrect & idxR & idxStimL & idxFR);
%     idxRed = find(idxC & idxCorrect & idxL & idxStimL & idxFR);
%     
%     thLeft = prctile(thMatrix(:, idxRed), prc, 2);
%     thRight = prctile(thMatrix(:, idxBlue), prc, 2);
%     plot(thLeft(:, 2), zAxis, 'r', 'LineWidth', 2)
%     hold on;
%     plot(thLeft(:, [1, 3]), zAxis, 'r:', 'LineWidth', 1);
%     plot(thRight(:, 2), zAxis, 'b', 'LineWidth', 2);
%     plot(thRight(:, [1, 3]), zAxis, 'b:', 'LineWidth', 1);
%     axis equal tight
%     ylim(zLims);
%     title(sprintf('Stim L (%1.0f; %1.0f)', length(idxRed), length(idxBlue)));
%     ylabel(sprintf('\\pm%1.0f %%', cc(iC)));
%     
%     subplot(nRows, 4, 2+(iC-2)*4)
%     idxBlue = find(idxC & idxCorrect & idxR & idxStimR & idxFR);
%     idxRed = find(idxC & idxCorrect & idxL & idxStimR & idxFR);
%     
%     thLeft = prctile(thMatrix(:, idxRed), prc, 2);
%     thRight = prctile(thMatrix(:, idxBlue), prc, 2);
%     plot(thLeft(:, 2), zAxis, 'r', 'LineWidth', 2)
%     hold on;
%     plot(thLeft(:, [1, 3]), zAxis, 'r:', 'LineWidth', 1);
%     plot(thRight(:, 2), zAxis, 'b', 'LineWidth', 2);
%     plot(thRight(:, [1, 3]), zAxis, 'b:', 'LineWidth', 1);
%     axis equal tight
%     ylim(zLims);
%     title(sprintf('Stim R (%1.0f; %1.0f)', length(idxRed), length(idxBlue)));
%     
%     subplot(nRows, 4, 3+(iC-2)*4)
%     idxBlue = find(idxC & idxCorrect & idxR & idxStimBoth & idxFR);
%     idxRed = find(idxC & idxCorrect & idxL & idxStimBoth & idxFR);
%     
%     thLeft = prctile(thMatrix(:, idxRed), prc, 2);
%     thRight = prctile(thMatrix(:, idxBlue), prc, 2);
%     plot(thLeft(:, 2), zAxis, 'r', 'LineWidth', 2)
%     hold on;
%     plot(thLeft(:, [1, 3]), zAxis, 'r:', 'LineWidth', 1);
%     plot(thRight(:, 2), zAxis, 'b', 'LineWidth', 2);
%     plot(thRight(:, [1, 3]), zAxis, 'b:', 'LineWidth', 1);
%     axis equal tight
%     ylim(zLims);
%     title(sprintf('Stim Both (%1.0f; %1.0f)', length(idxRed), length(idxBlue)));
%     
%     subplot(nRows, 4, 4+(iC-2)*4)
%     idxBlue = find(idxC & idxCorrect & idxR & idxStimNone & idxFR);
%     idxRed = find(idxC & idxCorrect & idxL & idxStimNone & idxFR);
%     
%     thLeft = prctile(thMatrix(:, idxRed), prc, 2);
%     thRight = prctile(thMatrix(:, idxBlue), prc, 2);
%     plot(thLeft(:, 2), zAxis, 'r', 'LineWidth', 2)
%     hold on;
%     plot(thLeft(:, [1, 3]), zAxis, 'r:', 'LineWidth', 1);
%     plot(thRight(:, 2), zAxis, 'b', 'LineWidth', 2);
%     plot(thRight(:, [1, 3]), zAxis, 'b:', 'LineWidth', 1);
%     axis equal tight
%     ylim(zLims);
%     title(sprintf('Stim None (%1.0f; %1.0f)', length(idxRed), length(idxBlue)));
%     
% end
% 
% %% plotting the trajectories - all in one
% 
% prc = [25 50 75];
% zLims = [5 60];
% cc = unique(abs(contrast));
% idxFR = finished & random;
% idxCorrect = outcome == 'C';
% idxR = behavior == 'R';
% idxL = behavior == 'L';
% idxStimL = optiStim(:, 1) > 0 & optiStim(:, 2) == 0;
% idxStimR = optiStim(:, 2) > 0 & optiStim(:, 1) == 0;
% idxStimBoth = optiStim(:, 1) > 0 & optiStim(:, 2) > 0;
% idxStimNone = optiStim(:, 1) == 0 & optiStim(:, 2) == 0;
% [zAxis, thMatrix] = binZTh(allZ, allTh, 30);
% 
% figure('Name', list(1).animalName)
% 
% nRows = 1;
% nC = length(cc);
% shade = [nC-1:-1:0]'/nC;
% for iC = 1:length(cc)
%     idxC = (abs(contrast) == cc(iC));
%     
%     subplot(nRows, 4, 1)
%     if cc(iC) == 0
%         idxBlue = find(idxC & idxR & idxStimL & idxFR);
%         idxRed = find(idxC & idxL & idxStimL & idxFR);
%     else
%         idxBlue = find(idxC & idxCorrect & idxR & idxStimL & idxFR);
%         idxRed = find(idxC & idxCorrect & idxL & idxStimL & idxFR);
%     end
%     
%     thLeft = prctile(thMatrix(:, idxRed), prc, 2);
%     thRight = prctile(thMatrix(:, idxBlue), prc, 2);
%     plot(thLeft(:, 2), zAxis, 'r', 'LineWidth', 2)
%     hold on;
%     plot(thLeft(:, [1, 3]), zAxis, 'r:', 'LineWidth', 1);
%     plot(thRight(:, 2), zAxis, 'b', 'LineWidth', 2);
%     plot(thRight(:, [1, 3]), zAxis, 'b:', 'LineWidth', 1);
%     axis equal tight
%     ylim(zLims);
%     title(sprintf('Stim L (%1.0f; %1.0f)', length(idxRed), length(idxBlue)));
%     ylabel(sprintf('\\pm%1.0f %%', cc(iC)));
%     
%     subplot(nRows, 4, 2+(iC-2)*4)
%     idxBlue = find(idxC & idxCorrect & idxR & idxStimR & idxFR);
%     idxRed = find(idxC & idxCorrect & idxL & idxStimR & idxFR);
%     
%     thLeft = prctile(thMatrix(:, idxRed), prc, 2);
%     thRight = prctile(thMatrix(:, idxBlue), prc, 2);
%     plot(thLeft(:, 2), zAxis, 'r', 'LineWidth', 2)
%     hold on;
%     plot(thLeft(:, [1, 3]), zAxis, 'r:', 'LineWidth', 1);
%     plot(thRight(:, 2), zAxis, 'b', 'LineWidth', 2);
%     plot(thRight(:, [1, 3]), zAxis, 'b:', 'LineWidth', 1);
%     axis equal tight
%     ylim(zLims);
%     title(sprintf('Stim R (%1.0f; %1.0f)', length(idxRed), length(idxBlue)));
%     
%     subplot(nRows, 4, 3+(iC-2)*4)
%     idxBlue = find(idxC & idxCorrect & idxR & idxStimBoth & idxFR);
%     idxRed = find(idxC & idxCorrect & idxL & idxStimBoth & idxFR);
%     
%     thLeft = prctile(thMatrix(:, idxRed), prc, 2);
%     thRight = prctile(thMatrix(:, idxBlue), prc, 2);
%     plot(thLeft(:, 2), zAxis, 'r', 'LineWidth', 2)
%     hold on;
%     plot(thLeft(:, [1, 3]), zAxis, 'r:', 'LineWidth', 1);
%     plot(thRight(:, 2), zAxis, 'b', 'LineWidth', 2);
%     plot(thRight(:, [1, 3]), zAxis, 'b:', 'LineWidth', 1);
%     axis equal tight
%     ylim(zLims);
%     title(sprintf('Stim Both (%1.0f; %1.0f)', length(idxRed), length(idxBlue)));
%     
%     subplot(nRows, 4, 4+(iC-2)*4)
%     idxBlue = find(idxC & idxCorrect & idxR & idxStimNone & idxFR);
%     idxRed = find(idxC & idxCorrect & idxL & idxStimNone & idxFR);
%     
%     thLeft = prctile(thMatrix(:, idxRed), prc, 2);
%     thRight = prctile(thMatrix(:, idxBlue), prc, 2);
%     plot(thLeft(:, 2), zAxis, 'r', 'LineWidth', 2)
%     hold on;
%     plot(thLeft(:, [1, 3]), zAxis, 'r:', 'LineWidth', 1);
%     plot(thRight(:, 2), zAxis, 'b', 'LineWidth', 2);
%     plot(thRight(:, [1, 3]), zAxis, 'b:', 'LineWidth', 1);
%     axis equal tight
%     ylim(zLims);
%     title(sprintf('Stim None (%1.0f; %1.0f)', length(idxRed), length(idxBlue)));
%     
% end
% 
% %%
% 
% cc = unique(abs(contrast));
% idxF = finished;
% idxRan = random;
% idxCorrect = outcome == 'C';
% idxR = behavior == 'R';
% idxL = behavior == 'L';
% idxStimL = optiStim(:, 1) > 0 & optiStim(:, 2) == 0;
% idxStimR = optiStim(:, 2) > 0 & optiStim(:, 1) == 0;
% idxStimBoth = optiStim(:, 1) > 0 & optiStim(:, 2) > 0;
% idxStimNone = optiStim(:, 1) == 0 & optiStim(:, 2) == 0;
% trialDur = nan(size(allZ));
% for iTrial = 1:length(allZ)
%     trialDur(iTrial) = sum(allZ{iTrial}>5)/60;
% end

% %%
% figure('Name', list(1).animalName)
% nFinBoth = sum(idxRan & idxF & idxStimBoth);
% nBoth = sum(idxRan & idxStimBoth);
% [pBoth, pciBoth] = binofit(nFinBoth, nBoth, alpha);
% nFinNone = sum(idxRan & idxF & idxStimNone); 
% nNone = sum(idxRan & idxStimNone);
% [pNone, pciNone] = binofit(nFinNone, nNone, alpha);
% nFinRight = sum(idxRan & idxF & idxStimR); 
% nRight = sum(idxRan & idxStimR);
% [pRight, pciRight] = binofit(nFinRight, nRight, alpha);
% nFinLeft = sum(idxRan & idxF & idxStimL); 
% nLeft = sum(idxRan & idxStimL);
% [pLeft, pciLeft] = binofit(nFinLeft, nLeft, alpha);
% 
% % histogram(trialDur(idxRan & idxStimNone))
% probs = [pNone, pRight, pLeft, pBoth];
% pciLow = probs - [pciNone(1), pciRight(1), pciLeft(1), pciBoth(1)];
% pciHigh = [pciNone(2), pciRight(2), pciLeft(2), pciBoth(2)] - probs;
% hBar = errorbar([1:4], probs, pciLow, pciHigh, 'o');
% ax = hBar.Parent;
% ax.XTick = 1:4;
% ax.XTickLabel = {'none', 'right', 'left', 'both'};
% xlim([0.5, 4.5]);
% title('Probability of finishing a random trial');
% xlabel('Inactivation condition');
% ylabel('Prob to finish with 95% c.i.');
% box off
% 
% figure('Name', list(1).animalName)
% meanBoth = mean(trialDur(idxRan & idxF & idxStimBoth));
% stdBoth = std(trialDur(idxRan & idxF & idxStimBoth))/sqrt(nFinBoth);
% meanNone = mean(trialDur(idxRan & idxF & idxStimNone));
% stdNone = std(trialDur(idxRan & idxF & idxStimNone))/sqrt(nFinNone);
% meanRight = mean(trialDur(idxRan & idxF & idxStimR));
% stdRight = std(trialDur(idxRan & idxF & idxStimR))/sqrt(nFinRight);
% meanLeft = mean(trialDur(idxRan & idxF & idxStimL));
% stdLeft = std(trialDur(idxRan & idxF & idxStimL))/sqrt(nFinLeft);
% 
% durs = [meanNone, meanRight, meanLeft, meanBoth];
% errs = [stdNone, stdRight, stdLeft, stdBoth];
% hBar = errorbar([1:4], durs, errs, 'o');
% ax = hBar.Parent;
% ax.XTick = 1:4;
% ax.XTickLabel = {'none', 'right', 'left', 'both'};
% xlim([0.5, 4.5]);
% title('Mean trial duration');
% xlabel('Inactivation condition');
% ylabel('trial duration [s] \pm s.e.m.');
% box off
% 
% figure('Name', list(1).animalName)
% [prcBoth] = prctile(trialDur(idxRan & idxF & idxStimBoth), [25 50 75]);
% [prcNone] = prctile(trialDur(idxRan & idxF & idxStimNone), [25 50 75]);
% [prcRight] = prctile(trialDur(idxRan & idxF & idxStimR), [25 50 75]);
% [prcLeft] = prctile(trialDur(idxRan & idxF & idxStimL), [25 50 75]);
% 
% durs = [prcNone(2), prcRight(2), prcLeft(2), prcBoth(2)];
% durLow = durs - [prcNone(1), prcRight(1), prcLeft(1), prcBoth(1)];
% durHigh = [prcNone(1), prcRight(1), prcLeft(1), prcBoth(1)] - durs;
% hBar = errorbar([1:4], durs, durLow, durHigh, 'o');
% ax = hBar.Parent;
% ax.XTick = 1:4;
% ax.XTickLabel = {'none', 'right', 'left', 'both'};
% xlim([0.5, 4.5]);
% title('Median \pm quartiles trial duration');
% xlabel('Inactivation condition');
% ylabel('trial duration [s]');
% box off
