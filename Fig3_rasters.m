ExpRef = '2014-08-15_1931_MK012';

figure('Name', ExpRef)

iPlane = 4;
iROI = 4;

%% plotting overfit rasters
nRows = 2;
nColumns = 3;

subplot(nRows, nColumns, 1)
plot([TM.modelEV.ZTh], [TM.modelEV.ZD], '.');
hold on;
plot(TM.modelEV(iPlane).ZTh(iROI), TM.modelEV(iPlane).ZD(iROI), 'or');
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('E.V. %')
xlabel('f(z, \theta)')
ylabel('f(z, d)')

subplot(nRows, nColumns, 2)
plot([TM.modelEV.ZTh], [TM.modelEV.ZThD], '.');
hold on;
plot(TM.modelEV(iPlane).ZTh(iROI), TM.modelEV(iPlane).ZThD(iROI), 'or');
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('E.V. %')
xlabel('f(z, \theta)')
ylabel('f(z, \theta, d)')

subplot(nRows, nColumns, 3)
plot([TM.modelEV.ZThD], [TM.modelEV.ZD], '.');
hold on;
plot(TM.modelEV(iPlane).ZThD(iROI), TM.modelEV(iPlane).ZD(iROI), 'or');
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('E.V. %')
xlabel('f(z, \theta, d)')
ylabel('f(z, d)')

subplot(nRows, nColumns, 4)
plot([TM.modelRho.ZTh], [TM.modelRho.ZD], '.');
hold on;
plot(TM.modelRho(iPlane).ZTh(iROI), TM.modelRho(iPlane).ZD(iROI), 'or');
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('\rho')
xlabel('f(z, \theta)')
ylabel('f(z, d)')

subplot(nRows, nColumns, 5)
plot([TM.modelRho.ZTh], [TM.modelRho.ZThD], '.');
hold on;
plot(TM.modelRho(iPlane).ZTh(iROI), TM.modelRho(iPlane).ZThD(iROI), 'or');
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('\rho')
xlabel('f(z, \theta)')
ylabel('f(z, \theta, d)')

subplot(nRows, nColumns, 6)
plot([TM.modelRho.ZThD], [TM.modelRho.ZD], '.');
hold on;
plot(TM.modelRho(iPlane).ZThD(iROI), TM.modelRho(iPlane).ZD(iROI), 'or');
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('\rho')
xlabel('f(z, \theta, d)')
ylabel('f(z, d)')

figure('Name', ExpRef)
subplot(1, 3, 1)
plot([TM.modelEV.ZTh], [TM.modelRho.ZTh], '.')
hold on;
plot(TM.modelEV(iPlane).ZTh(iROI), TM.modelRho(iPlane).ZTh(iROI), 'or');
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('f(z, \theta)')
xlabel('E.V.')
ylabel('\rho')

subplot(1, 3, 2)
plot([TM.modelEV.ZD], [TM.modelRho.ZD], '.')
hold on;
plot(TM.modelEV(iPlane).ZD(iROI), TM.modelRho(iPlane).ZD(iROI), 'or');
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('f(z, d)')
xlabel('E.V.')
ylabel('\rho')

subplot(1, 3, 3)
plot([TM.modelEV.ZThD], [TM.modelRho.ZThD], '.')
hold on;
plot(TM.modelEV(iPlane).ZThD(iROI), TM.modelRho(iPlane).ZThD(iROI), 'or');
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('f(z, \theta, d)')
xlabel('E.V.')
ylabel('\rho')

%% plotting cross-validated rasters
figure('Name', ExpRef)

nRows = 2;
nColumns = 3;

subplot(nRows, nColumns, 1)
plot([TM.modelEV.ZThCV], [TM.modelEV.ZDCV], '.');
hold on;
plot(TM.modelEV(iPlane).ZThCV(iROI), TM.modelEV(iPlane).ZDCV(iROI), 'or');
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('E.V. %')
xlabel('f(z, \theta)')
ylabel('f(z, d)')
box off;

subplot(nRows, nColumns, 2)
plot([TM.modelEV.ZThCV], [TM.modelEV.ZThDCV], '.');
hold on;
plot(TM.modelEV(iPlane).ZThCV(iROI), TM.modelEV(iPlane).ZThDCV(iROI), 'or');
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('E.V. %')
xlabel('f(z, \theta)')
ylabel('f(z, \theta, d)')
box off;

subplot(nRows, nColumns, 3)
plot([TM.modelEV.ZThDCV], [TM.modelEV.ZDCV], '.');
hold on;
plot(TM.modelEV(iPlane).ZThDCV(iROI), TM.modelEV(iPlane).ZDCV(iROI), 'or');
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('E.V. %')
xlabel('f(z, \theta, d)')
ylabel('f(z, d)')
box off;

subplot(nRows, nColumns, 4)
plot([TM.modelRho.ZThCV], [TM.modelRho.ZDCV], '.');
hold on;
plot(TM.modelRho(iPlane).ZThCV(iROI), TM.modelRho(iPlane).ZDCV(iROI), 'or');
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('\rho')
xlabel('f(z, \theta)')
ylabel('f(z, d)')
box off;

subplot(nRows, nColumns, 5)
plot([TM.modelRho.ZThCV], [TM.modelRho.ZThDCV], '.');
hold on;
plot(TM.modelRho(iPlane).ZThCV(iROI), TM.modelRho(iPlane).ZThDCV(iROI), 'or');
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('\rho')
xlabel('f(z, \theta)')
ylabel('f(z, \theta, d)')
box off;

subplot(nRows, nColumns, 6)
plot([TM.modelRho.ZThDCV], [TM.modelRho.ZDCV], '.');
hold on;
plot(TM.modelRho(iPlane).ZThDCV(iROI), TM.modelRho(iPlane).ZDCV(iROI), 'or');
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('\rho')
xlabel('f(z, \theta, d)')
ylabel('f(z, d)')
box off;

figure('Name', ExpRef)
subplot(1, 3, 1)
plot([TM.modelEV.ZThCV], [TM.modelRho.ZThCV], '.')
hold on;
plot(TM.modelEV(iPlane).ZThCV(iROI), TM.modelRho(iPlane).ZThCV(iROI), 'or');
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('f(z, \theta)')
xlabel('E.V.')
ylabel('\rho')

subplot(1, 3, 2)
plot([TM.modelEV.ZDCV], [TM.modelRho.ZDCV], '.')
hold on;
plot(TM.modelEV(iPlane).ZDCV(iROI), TM.modelRho(iPlane).ZDCV(iROI), 'or');
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('f(z, d)')
xlabel('E.V.')
ylabel('\rho')

subplot(1, 3, 3)
plot([TM.modelEV.ZThDCV], [TM.modelRho.ZThDCV], '.')
hold on;
plot(TM.modelEV(iPlane).ZThDCV(iROI), TM.modelRho(iPlane).ZThDCV(iROI), 'or');
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('f(z, \theta, d)')
xlabel('E.V.')
ylabel('\rho')

%% comparison of overfitted and cross-validated E.V. and rho

nRows = 2;
nColumns = 3;

subplot(nRows, nColumns, 1)
plot([TM.modelEV.ZTh], [TM.modelEV.ZThCV], '.');
hold on;
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('f(z, \theta) E.V. %')
xlabel('Overfit')
ylabel('Cross-validated')

subplot(nRows, nColumns, 2)
plot([TM.modelEV.ZThD], [TM.modelEV.ZThDCV], '.');
hold on;
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('f(z, \theta, d) E.V. %')
xlabel('Overfit')
ylabel('Cross-validated')

subplot(nRows, nColumns, 3)
plot([TM.modelEV.ZD], [TM.modelEV.ZDCV], '.');
hold on;
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('f(z, d) E.V. %')
xlabel('Overfit')
ylabel('Cross-validated')

subplot(nRows, nColumns, 4)
plot([TM.modelRho.ZTh], [TM.modelRho.ZThCV], '.');
hold on;
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('f(z, \theta) \rho')
xlabel('Overfit')
ylabel('Cross-validated')

subplot(nRows, nColumns, 5)
plot([TM.modelRho.ZThD], [TM.modelRho.ZThDCV], '.');
hold on;
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('f(z, \theta, d) \rho')
xlabel('Overfit')
ylabel('Cross-validated')

subplot(nRows, nColumns, 6)
plot([TM.modelRho.ZD], [TM.modelRho.ZDCV], '.');
hold on;
plot([0 1], [0 1], 'k:');
axis equal square;
xlim([0 1]);
ylim([0 1])
title('f(z, d) \rho')
xlabel('Overfit')
ylabel('Cross-validated')
