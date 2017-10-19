function nPlanes=getNumberOfPlanes(info)

% this function analizes the Timeline data and extracts the actual number
% of planes in the dataset.
% it uses the Piezo Command signal and the Frame Trigger signal
% Currently can be used only with 'Sawtooth' pieo command waveform, which
% goes from low to high values (from shallow to deep layers, and then jumps
% back to shallow)

% 2014-02-24 - MK created

% try
%     load(fullfile(info.folderTLLocal, [info.basenameTL, '.mat']));
% catch
%     fprintf('Loading Timeline data from zserver2 (local loading failed)... \n');
%     load(fullfile(info.folderTL, [info.basenameTL, '.mat']));
% end
% 
% nInputs=length(Timeline.hw.inputs);
% indFrames=nan;
% indPiezo=nan;
% for iInput=1:nInputs
%     if isequal(Timeline.hw.inputs(iInput).name, 'neuralFrames')
%         indFrames=iInput;
%     elseif isequal(Timeline.hw.inputs(iInput).name, 'piezoCommand')
%         indPiezo=iInput;
%     else
%     end
% end
% 
% framePulses=diff(Timeline.rawDAQData(:, indFrames));
% maxlag = min(round(length(framePulses)/10), Timeline.hw.daqSampleRate*300);
% xFramePulses = xcorr(framePulses, maxlag, 'unbiased');
% xFramePulses = filtfilt([1 1 1]/3, 1, xFramePulses);
% [pks, locs] = findpeaks(xFramePulses);
% th = 0.1*max(pks);
% pkLocs = locs(pks>th);
% frameRate = 1/mean(diff(pkLocs))*Timeline.hw.daqSampleRate;
% 
% % % plotting for debugging
% % figure;
% % subplot(2, 1, 1);
% % plot(xFramePulses);
% % hold on;
% % plot(locs, pks, 'r.', xlim, [th th], 'c:');
% 
% piezoSignal=diff(Timeline.rawDAQData(:, indPiezo));
% maxlag = min(round(length(piezoSignal)/10), Timeline.hw.daqSampleRate*300);
% xPiezoSignal = xcorr(piezoSignal, maxlag, 'unbiased');
% xPiezoSignal = filtfilt([1 1 1]/3, 1, xPiezoSignal);
% [pks, locs] = findpeaks(xPiezoSignal);
% th = 0.2*max(pks);
% pkLocs = locs(pks>th);
% cycleRate = 1/mean(diff(pkLocs))*Timeline.hw.daqSampleRate;
% 
% % % plotting for debugging
% % subplot(2, 1, 2);
% % plot(xPiezoSignal);
% % hold on;
% % plot(locs, pks, 'r.', xlim, [th th], 'c:');
% 
% nPlanesEst = frameRate/cycleRate;
% % 1-mod(nPlanes,1)
% 
% if nPlanesEst~=round(nPlanesEst)
%     fprintf('Estimated nPlanes = %12.10f\n', nPlanesEst);
%     warning('Uncertain about the estimate from the piezo signal, number of planes is not an integer');
%     fprintf('Will now get the number of planes from the tiff header...\n');
% end

try
    allTiffInfo = dir([info.folder2pLocal, filesep, info.basename2p, '*.tif']);
    tiffName = allTiffInfo(1).name;
    filename=fullfile(info.folder2pLocal, tiffName);
    [~, header]=img.loadFrames(filename, 1, 1, 1);
catch
    fprintf('Getting the tiff from the server (local tiffs do not exist)...\n');
    allTiffInfo = dir([info.folder2p, filesep, info.basename2p, '*.tif']);
    tiffName = allTiffInfo(1).name;
    filename=fullfile(info.folder2p, tiffName);
    [~, header]=img.loadFrames(filename, 1, 1, 1);
end
% getting some parameters from the header
hh=header{1};
fastZEnable = sscanf(hh(findstr(hh, 'fastZEnable = '):end), 'fastZEnable = %d');
fastZDiscardFlybackFrames = sscanf(hh(findstr(hh, 'fastZDiscardFlybackFrames = '):end), 'fastZDiscardFlybackFrames = %d');
if isempty(fastZDiscardFlybackFrames)
    fastZDiscardFlybackFrames = 0;
end
stackNumSlices = sscanf(hh(findstr(hh, 'stackNumSlices = '):end), 'stackNumSlices = %d');

if fastZEnable
    nPlanes=stackNumSlices+fastZDiscardFlybackFrames
%     if (nPlanesEst-nPlanes)/nPlanes>1e-4
%         warning('The difference between the piezo-signal estimated and the ScanImage number of planes is large. Check your data');
%     end
else
    fprintf('The fast scanning was disabled during this acquisition\n');
    nPlanes=1
end;




