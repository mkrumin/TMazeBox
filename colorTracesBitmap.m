load('G:\DATA\2014-08-15_1931_MK012_TM.mat')
% load('G:\DATA\2017-09-23_1539_JL008_TM.mat')
% load('G:\DATA\2017-09-24_1558_JL008_TM.mat')


res = TM.trainingData;
nPlanes = length(res);
ev = [];
planeInd = [];
roiInd = [];
for iPlane = 1:nPlanes
    nROIs = length(res{iPlane});
    planeInd = cat(1, planeInd, ones(nROIs, 1)*iPlane);
    roiInd = cat(1, roiInd, [1:nROIs]');
    ev = cat(1, ev, 1-[res{iPlane}.errVals]');
end

[~, sortedIdx] = sort(ev, 'descend');

for iCell = 1:1
%     iPlane = planeInd(sortedIdx(iCell));
%     iROI = roiInd(sortedIdx(iCell));
    
    iPlane = 4;
    iROI = 4;
    
    map = res{iPlane}(iROI).zThetaMap;
    zAxis = res{iPlane}(iROI).zThetaBinCentres{1};
    thAxis = res{iPlane}(iROI).zThetaBinCentres{2};
    
    fData = TM.data2p{iPlane}.F(:, iROI);
    f0 = prctile(fData, 20);
    fData = (fData-f0)/f0;
    tData = TM.times2p{iPlane}';
    
    contrasts = TM.dataTMaze.contrastSequence;
    report = TM.dataTMaze.report';
    cc = unique(contrasts);
    nRows = 2;
    nColumns = length(cc);
    
    nTrials = TM.dataTMaze.nTrials;
    hFig = figure(1);
    clf;
    hFig.Color = [1 1 1];
    
    for iTrial = 1:nTrials
        if ~ismember(report(iTrial), 'RL')
            % we only want finished trials
            continue;
        end
        t = TM.timesVRframes(iTrial).t(2:end-1);
        if isempty(t)
            continue;
        end
        idx = TM.timesVRframes(iTrial).idx(2:end);
        z = -TM.dataTMaze.SESSION.allTrials(iTrial).posdata(idx,3);
        theta = 180/pi*TM.dataTMaze.SESSION.allTrials(iTrial).posdata(idx,4);
        f = interp1(tData, fData, t);
        
        iRow = 1*(report(iTrial) == 'L') + 2*(report(iTrial) == 'R');
        iColumn = find(cc == contrasts(iTrial));
        iPlot = sub2ind([nColumns, nRows], iColumn, iRow);
        side(iPlot) = report(iTrial);
        contr(iPlot) = contrasts(iTrial);
        %     if iPlot ~= 1
        %         continue;
        %     end
        
        ax(iPlot) = subplot(nRows, nColumns, iPlot);
        theta = theta(1:10:end);
        z = z(1:10:end);
        f = f(1:10:end);
        %     s = surface([theta, theta]', [z, z]', 0*[z, z]', [f, f]');
        %     s.EdgeColor = 'interp';
        %     s.LineWidth = 1;
        
        theta(end+1) = NaN;
        z(end+1) = NaN;
        f(end+1) = NaN;
        p = patch(theta, z, f);
        p.EdgeColor = 'interp';
        p.LineWidth = 2;
        
        hold on;
        if iRow == 1
            title(sprintf('%1.0f%%', contrasts(iTrial)));
            set(gca, 'XTickLabel', '')
        end
        if iRow == 2 && iColumn == 1
            lbl = xlabel('\theta [deg]');
            lbl.FontSize = 12;
        end
        if iColumn == 1
            lbl = ylabel('z [cm]');
            lbl.FontSize = 12;
        else
            set(gca, 'YTickLabel', '')
            set(gca, 'XTickLabel', '')
        end
    end
    
    cminmax = prctile(fData, [0.1 99.9]);
%     cminmax = [0 2]; % to match the existing figure
%     thr = mean(prctile(fData, [0.1 99.9]));
    % thr = 1;
    thr = min(map(:)) + 0.6*(max(map(:))-min(map(:)));
    cont = contourc(thAxis, zAxis, map, [thr thr]);
    cX = cont(1, 2:end);
    cX(end+1) = cX(1);
    cY = cont(2, 2:end);
    cY(end+1) = cY(1);
    
%     ax = hFig.Children;
    for iAx = 1:length(ax)
        colormap(ax(iAx), 'jet');
        caxis(ax(iAx), cminmax);
        %     ax(iAx).XLim = ([min(thAxis), max(thAxis)]);
        ax(iAx).XLim = ([-30, 30]);
        ax(iAx).YLim = ([0, max(zAxis)]);
        ax(iAx).DataAspectRatio = [1 1 1];
        ax(iAx).XTick = [-30 0 30];
        ax(iAx).YTick = [0 50 100];
        ax(iAx).FontSize = 14;
        %     plot(ax(iAx), cX, cY, 'k--', 'LineWidth', 2);
    end
    
%     pause;
end

%%

h = figure('Position', [200 200 344 600]);
h.Color = [1 1 1];
for iAx = 1:length(ax)
a = ax(iAx);
oldParent = a.Parent;
oldPosition = a.Position;
a.Position = oldPosition; % to make it 'manual'?
a.Parent = h;
pause(1);
a.Position = [0 0 1 1];
pause(1);
a.Visible = 'off';
pause(1);
filename = sprintf('%s_%2.0f.png', side(iAx), contr(iAx));
print(h, fullfile('G:\Temp2p\', filename), '-dpng');
% filename = sprintf('%s_%2.0f.jpg', side(iAx), contr(iAx));
% print(h, fullfile('G:\Temp2p\', filename), '-djpeg');
% fr = getframe(a);
% cla(a);
a.Parent = oldParent;
a.Position = oldPosition;
% [nY, nX, ~]  = size(fr.cdata);
% xl = a.XLim;
% xAxis = linspace(xl(1), xl(2), nX);
% yl = a.YLim;
% yAxis = linspace(yl(1), yl(2), nY);
% imagesc(a, xAxis, yAxis, flipud(fr.cdata));
% a.YDir = 'normal';
a.Visible = 'on';
end
