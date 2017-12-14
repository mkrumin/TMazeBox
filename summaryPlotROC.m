warning off
folder = 'G:\DATA\';

files = dir('G:\DATA\*_TM.mat');
nFiles = length(files);
nRows = floor(sqrt(nFiles));
nColumns = ceil(nFiles/nRows);
sigFrac = zeros(nFiles, 2);
figure
for iFile = 1:nFiles
    load(fullfile(folder, files(iFile).name));
    subplot(nRows, nColumns, iFile);
    sigFrac(iFile, :) = TM.rocAnalysisTmp;
    drawnow;
end

figure
plot(sigFrac(:,1), sigFrac(:,2), 'o');
hold on;
plot([0, 1], [0, 1], 'k:');
xlabel('Full Data');
ylabel('Residuals');
title('Fraction of cells with significant AUROC');
warning on
