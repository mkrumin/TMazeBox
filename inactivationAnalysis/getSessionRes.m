function out = getSessionRes(EXP, SESSION)

nTrials = length(SESSION.allTrials);
try
    out.ExpRef = EXP.expRef;
catch
    out.ExpRef = '';
end
out.contrast = nan(nTrials, 1);
out.outcome = '';
out.behavior = '';
out.finished = nan(nTrials, 1);
out.random = nan(nTrials, 1);
% HACK assuming two locations
out.optiStim = nan(nTrials, 2);
out.isV1 = nan(nTrials, 1);
out.isPPC = nan(nTrials, 1);
out.isUndecided = nan(nTrials, 1);
out.z = cell(nTrials, 1);
out.theta  = cell(nTrials, 1);

for iTrial = 1:nTrials
    out.contrast(iTrial) = SESSION.allTrials(iTrial).info.contrast;
    side = -1+2*isequal(SESSION.allTrials(iTrial).info.stimulus, 'RIGHT');
    out.contrast(iTrial) = out.contrast(iTrial) * side;
    out.outcome(iTrial) = SESSION.allTrials(iTrial).info.outcome(1);
    out.behavior(iTrial) = out.outcome(iTrial);
    if out.outcome(iTrial) == 'C'
        out.behavior(iTrial) = SESSION.allTrials(iTrial).info.stimulus(1);
    elseif out.outcome(iTrial) == 'W'
        out.behavior(iTrial) = char('R'+'L'-SESSION.allTrials(iTrial).info.stimulus(1));
    end
    out.finished(iTrial) = ismember(out.behavior(iTrial), {'R', 'L'});
    if isequal(EXP.stimType, 'BAITED')
        out.random(iTrial) = iTrial == 1 || out.outcome(iTrial-1)=='C';
    else
        out.random(iTrial) = true; %assuming 'RANDOM'
    end
    out.optiStim(iTrial, :) = SESSION.allTrials(iTrial).info.optiStim.laserPower;
    
    
    out.isV1(iTrial) = isequal(SESSION.allTrials(iTrial).info.optiStim.ML, [-2.5, 2.5]);
    out.isPPC(iTrial) = isequal(SESSION.allTrials(iTrial).info.optiStim.ML, [-1.7, 1.7]);
    
    zInd = find(ismember(SESSION.allTrials(iTrial).pospars, 'Z'));
    thInd = find(ismember(SESSION.allTrials(iTrial).pospars, 'theta'));
    out.z{iTrial} = -SESSION.allTrials(iTrial).posdata(:, zInd);
    out.theta{iTrial} = SESSION.allTrials(iTrial).posdata(:, thInd)*180/pi;
end
zOptiStimStart = max(10, EXP.optiStimList(1).onset);
[out.isUndecided, out.thPreStim] = isUndecidedPrestim(out, zOptiStimStart);
out.behavior = out.behavior(:);
out.outcome = out.outcome(:);
