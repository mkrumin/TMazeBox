function out = getTrialDur(TM)


nTrials = TM.dataTMaze.nTrials;
dur = nan(nTrials, 1);
for iTrial = 1:nTrials
    
    frameTimes = TM.dataTMaze.SESSION.allTrials(iTrial).time(2:end);
    idxStart = find(TM.dataTMaze.SESSION.allTrials(iTrial).freezeOver, 1, 'first');
    dur(iTrial) = frameTimes(end) - frameTimes(idxStart);
end

idxFinished = ismember(TM.dataTMaze.report, 'RL')';

% figure
% histogram(dur(idxFinished))

out.ExpRef = TM.expRef;
out.dur = dur;
out.isFinished = idxFinished;
