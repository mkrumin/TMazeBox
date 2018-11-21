% this script will generate datasets to be shared with the eLife paper
% Data to be included:
%         - Fluorescence traces after Suite2p
%         - x, z, theta traces
%         - traces of three components of the ball rotation
%         - identities of the trials - contrast, reward


folder = 'G:\DATA\';

allExpRefs = {...
    '2014-08-02_2023_MK012'; ...
    '2014-08-04_1827_MK012'; ...
    '2014-08-05_2228_MK012'; ...
    '2014-08-08_2251_MK012'; ...
    '2014-08-11_2133_MK012'; ...
    '2014-08-13_2023_MK012'; ...
    '2014-08-15_1931_MK012'; ...
    '2014-08-02_2203_MK014'; ...
    '2014-08-05_1937_MK014'; ...
    '2017-06-12_1420_JL005'; ...
    '2017-07-15_1708_JL008'; ...
    '2017-07-27_1433_JL008'; ...
    '2017-08-12_1056_JL008'; ...
    '2017-09-18_1707_JL008'; ...
    '2017-09-23_1539_JL008'; ...
    '2015-07-03_2127_MK020'; ...
    '2015-07-07_2127_MK020'; ...
    '2015-07-30_2010_MK020'; ...
    '2015-07-31_1851_MK020'; ...
    '2015-08-02_1709_MK020'; ...
    '2015-08-05_1930_MK020'; ...
    '2015-08-07_1821_MK020'; ...
    '2015-07-03_1632_MK022'; ...
    '2015-08-02_2051_MK022'; ...
    '2015-06-19_2155_MK023'; ...
    '2015-07-18_154_MK023';
    };

allAnimals = {'MK012', 'MK014', 'MK020', 'MK022', 'MK023', 'JL005', 'JL008'};
allGenotypes = [repmat({'C57bl/6'}, 1, 2), ...
    repmat({'Camk2a-tTA; Ai93(TITL-GCaMP6f); Emx1-IRES-Cre'}, 1, 3), ...
    repmat({'Ai95(RCL-GCaMP6f)-D; Slc17a7-IRES2-Cre-D'}, 1, 2)];

%%

allExpRefs = allExpRefs(15);
nExperiments = length(allExpRefs);
% nExperiments = 1;

data = struct('subject', '', 'genotype', '', 'nPlanes', [], 'traces', [], 'ball', [], 'trials', []);
for iExp = 1:nExperiments
    fprintf('Processing dataset %g/%g...\n', iExp, nExperiments);
    % create TM object
    TM = TMazeVR(allExpRefs{iExp});
    % or load already existing data from disk
%     load(fullfile(folder, [allExpRefs{iExp}, '_TM.mat']));
    ball.t = TM.dataBall.t;
    ball.pitch = TM.dataBall.forward;
    ball.roll = TM.dataBall.sideways;
    ball.yaw = -TM.dataBall.rotation;
    nTrials = TM.dataTMaze.nTrials;
    trials = struct( 'contrast', [], 'report', '', 'rewarded', [], 't', [], 'x', [], 'z', [], 'theta', []);
    pospars = TM.dataTMaze.SESSION.allTrials(1).pospars;
    [~, xInd] = ismember('X', pospars);
    [~, zInd] = ismember('Z', pospars);
    [~, thInd] = ismember('theta', pospars);
    
    for iTrial = 2:nTrials
        idx = TM.timesVRframes(iTrial).idx;
        idx = idx(2:end);
        trials(iTrial).t = TM.timesVRframes(iTrial).t(2:end-1);
        trials(iTrial).x = TM.dataTMaze.SESSION.allTrials(iTrial).posdata(idx, xInd);
        trials(iTrial).z = - TM.dataTMaze.SESSION.allTrials(iTrial).posdata(idx, zInd);
        trials(iTrial).theta = TM.dataTMaze.SESSION.allTrials(iTrial).posdata(idx, thInd);
        trials(iTrial).contrast = TM.dataTMaze.contrastSequence(iTrial);
        trials(iTrial).report = TM.dataTMaze.report(iTrial);
        trials(iTrial).rewarded = TM.dataTMaze.outcome(iTrial) == 'C';
    end
    
    %%
    fprintf('Estimating optimal scaling factors..');
    dt = 0.001;
    maxlag = 1000;
    scalingZ = [];
    scalingTh = [];
    rhoZ = [];
    rhoTh = [];
    
    for iTrial = 2:nTrials
        tStart = ceil(trials(iTrial).t(1)/dt)*dt+0.1;
        tEnd = floor(trials(iTrial).t(end)/dt)*dt-0.1;
        tt = tStart:dt:tEnd;
        
        run = interp1(ball.t, ball.pitch, tt);
        dz = diff(trials(iTrial).z) ./ cos(trials(iTrial).theta(2:end));
        dz = interp1(trials(iTrial).t(1:end-1), medfilt1(dz, 5), tt);
        
        yaw = interp1(ball.t, ball.yaw, tt);
        dth = diff(trials(iTrial).theta);
        dth = interp1(trials(iTrial).t(1:end-1), medfilt1(dth), tt);
        
        xcZ = xcorr(run, dz, maxlag);
        [mxcZ, mIndZ] = max(xcZ);
        deltaT = mIndZ-maxlag-1;
        if deltaT<=0
            dzShifted = dz(1 - deltaT : end);
            runShifted = run(1 : end + deltaT);
        elseif deltaT>0
            dzShifted = dz(1 : end - deltaT);
            runShifted = run(1 + deltaT : end);
        end
        ttShifted = tt(1:end-abs(deltaT));
        scalingZ(iTrial) = dzShifted(:)\runShifted(:);
        tmp = corrcoef(runShifted, dzShifted);
        rhoZ(iTrial) = tmp(2);
        
%         figure
%         subplot(2, 2, 1)
%         plot(-maxlag:maxlag, xcZ)
%         hold on;
%         plot([0 0], ylim, 'k:');
%         plot(mIndZ-maxlag-1, mxcZ, '*r');
%         text(mIndZ-maxlag-1, mxcZ, sprintf('\\Deltat = %01.0f ms', mIndZ-maxlag-1))
%         
%         ax1 = subplot(2, 2, 2);
%         plot(ttShifted, runShifted, ttShifted, dzShifted*scalingZ(iTrial));
%         title(rhoZ(iTrial))
        
        xcTh = xcorr(yaw, dth, maxlag);
        [mxcTh, mIndTh] = max(xcTh);
        deltaT = mIndTh-maxlag-1;
        if deltaT<=0
            dthShifted = dth(1 - deltaT : end);
            yawShifted = yaw(1 : end + deltaT);
        elseif deltaT>0
            dthShifted = dth(1 : end - deltaT);
            yawShifted = yaw(1 + deltaT : end);
        end
        ttShifted = tt(1:end-abs(deltaT));
        scalingTh(iTrial) = dthShifted(:)\yawShifted(:);
        tmp = corrcoef(yawShifted, dthShifted);
        rhoTh(iTrial) = tmp(2);
        
%         subplot(2, 2, 3)
%         plot(-maxlag:maxlag, xcTh)
%         hold on;
%         plot([0 0], ylim, 'k:');
%         plot(mIndTh-maxlag-1, mxcTh, '*r');
%         text(mIndTh-maxlag-1, mxcTh, sprintf('\\Deltat = %01.0f ms', mIndTh-maxlag-1))
%         
%         ax2 = subplot(2, 2, 4);
%         plot(ttShifted, yawShifted, ttShifted, dthShifted*scalingTh(iTrial));
%         title(rhoTh(iTrial))
%         
%         linkaxes([ax1, ax2], 'x')
    end
    fprintf('.done\n');
    scZ = mean(scalingZ(rhoZ>0.97));
    scTh = mean(scalingTh(rhoTh>0.97));
    
    dtBall = median(diff(ball.t));
    ball.pitch = ball.pitch / scZ / dtBall;
    ball.roll = ball.roll / scZ / dtBall;
    ball.yaw = ball.yaw / scTh / dtBall;
    %%
    %     keyboard;
    
    % testing the scaling and alignment by plotting
%     figure('Name', allExpRefs{iExp});
%     subplot(2, 1, 1)
%     plot(ball.t, cumsum(ball.yaw)*dtBall)
%     hold on;
%     for iTrial = 2:nTrials 
%         tmp = cumsum(ball.yaw)*dtBall; 
%         [~, ind] = min(abs(ball.t - trials(iTrial).t(1))); 
%         shift = tmp(ind); 
%         plot(trials(iTrial).t, trials(iTrial).theta + shift, 'LineWidth', 3); 
%     end
%     xlabel('time [sec]')
%     ylabel('Cumulative forward running [cm]');
%     
%     subplot(2, 1, 2)
%     plot(ball.t, cumsum(ball.pitch)*dtBall)
%     hold on;
%     for iTrial = 2:nTrials 
%         tmp = cumsum(ball.pitch)*dtBall; 
%         [~, ind] = min(abs(ball.t - trials(iTrial).t(1))); 
%         shift = tmp(ind); 
%         plot(trials(iTrial).t, trials(iTrial).z + shift, 'LineWidth', 3); 
%     end
%     xlabel('time [sec]')
%     ylabel('Cumulative heading angle [rad]');
    
    % now let's do the fluorescence data
    planes = TM.Planes;
    traces = struct('nCells', [], 't', [], 'F', []);
    for iPlane = 1:length(planes)
        validIdx = ismember(TM.data2p{planes(iPlane)}.ROI.CellClasses, 's');
        F = TM.data2p{planes(iPlane)}.F(:, validIdx);
        validIdx = ~all(isnan(F));
        F = F(:, validIdx);
        t = TM.times2p{planes(iPlane)}(:);
        traces(iPlane).F = F;
        traces(iPlane).t = t;
        traces(iPlane).nCells = size(traces(iPlane).F, 2);
    end
    
    % excluding the first (empty) trial
    trials = trials(2:end);

    data(iExp).subject = allExpRefs{iExp}(end-4:end);
    data(iExp).nPlanes = length(traces);
    [~, idx] = ismember(data(iExp).subject, allAnimals); 
    data(iExp).genotype = allGenotypes{idx};
    data(iExp).traces = traces(:);
    data(iExp).ball = ball;
    data(iExp).trials = trials(:);
end

data = data(:);

%%
% fprintf('Saving data to disk..');
% save(fullfile(folder, 'Krumin_etal_2018_eLife.mat'), 'data');
% fprintf('.done! YAY!\n');