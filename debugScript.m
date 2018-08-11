% debug script

figure, 

nRows = 2;
nColumns = 3;

% trajectories

subplot(nRows, nColumns, 1);
for iTrial = 1:nTrials 
    plot(th{iTrial}, zz{iTrial}); 
    hold on; 
end, 
axis equal

xlim([-30 30]);
ylim([5 105]);

% models

imagesc(thAxis, zAxis, oldMap); 
colorbar; 
axis equal tight xy;
