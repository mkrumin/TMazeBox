function hPlot = plotSbySSummary(hAx, xx, yy, plotType, axesLims)

switch plotType
    case 'mean'
        hPlot = plot(mean(xx), mean(yy), 'o');
    case 'median'
        hPlot = plot(median(xx), median(yy), 'o');
    case 'meanstd'
        xMean = mean(xx);
        xNeg = std(xx);
        xPos = xNeg;
        yMean = mean(yy);
        yNeg = std(yy);
        yPos = yNeg;
        
        er = errorbar(xData, yData, yNeg, yPos, xNeg, xPos, 'o');
        er.MarkerSize = 6;
        er.MarkerFaceColor = er.MarkerEdgeColor;
        er.LineStyle = ':';
        er.LineWidth = 0.25;
        hPlot = er;
    case 'medianprc'
        xPr = prctile(xx, [25 50 75]);
        xData = xPr(2);
        xNeg = xPr(2) - xPr(1);
        xPos = xPr(3) - xPr(2);
        yPr = prctile(yy, [25 50 75]);
        yData = yPr(2);
        yNeg = yPr(2) - yPr(1);
        yPos = yPr(3) - yPr(2);
        
        er = errorbar(xData, yData, yNeg, yPos, xNeg, xPos, 'o');
        er.MarkerSize = 6;
        er.MarkerFaceColor = er.MarkerEdgeColor;
        er.LineStyle = ':';
        er.LineWidth = 0.25;
        hPlot = er;
    otherwise
end

hold on;
plot(axesLims, axesLims, 'k--');
axis equal
xlim(axesLims)
ylim(axesLims)
% plot([0 0], axesLims, 'k:', axesLims, [0 0], 'k:');
