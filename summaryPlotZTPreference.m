warning off
folder = 'G:\DATA\';

files = dir('G:\DATA\*_TM.mat');
nFiles = length(files);
nRows = floor(sqrt(nFiles));
nColumns = ceil(nFiles/nRows);
sigFrac = zeros(nFiles, 2);
figure
for iFile = 1:nFiles
    fprintf('Analyzing file %d/%d\n', iFile, nFiles)
    load(fullfile(folder, files(iFile).name));
    subplot(nRows, nColumns, iFile);
    [sigFrac(iFile, :), auroc{iFile}] = TM.rocAnalysisTmp;
    drawnow;
end

return
%%
figure
plot(sigFrac(:,1), sigFrac(:,2), 'o');
hold on;
plot([0, 1], [0, 1], 'k:');
axis equal tight
box off
xlabel('Full Data');
ylabel('Residuals');
title('Fraction of cells with significant AUROC');
set(gca, 'XTick', [0 0.5 1], 'YTick', [0 0.5 1]);
warning on


%%
figure
aurocMat = cell2mat(auroc(:));
plot(aurocMat(:,1), aurocMat(:,2), '.');
xlim([0 1]);
ylim([0 1]);
axis equal 
set(gca, 'XTick', [0 0.5 1], 'YTick', [0 0.5 1]);
xlabel('Full Data AUROC');
ylabel('Residuals AUROC');
box off

figure;
histogram(aurocMat(:,1));
hold on;
histogram(aurocMat(:,2));
legend('Full data', 'Residuals');
xlabel('AUROC');
box off