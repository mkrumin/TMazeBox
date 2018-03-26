function getAllExps(subject)

if nargin<1
    subject = 'JC001';
    %     subject = 'M140711_MK014';
end
currDir = cd('\\zserver.cortexlab.net\Code\Rigging\main');
oldList = dat.listExps(subject);

for iExp = 1:length(oldList)
    
    try
        tm = TMaze(oldList{iExp});
        fprintf('%s: % 4.0f trials\n', oldList{iExp}, tm.nTrials);
    catch
    end
end

cd(currDir);
newList = dat.listExps(subject);
newList = newList(end-1:end);
nExps = length(newList);
for iExp = nExps-1:nExps
    
    try
        tm = TMaze(newList{iExp});
        fprintf('%s: % 4.0f trials\n', newList{iExp}, tm.nTrials);
        tm = tm.getPCData;
        tm = tm.fitPC;
        ax = subplot(1, nExps, iExp);
        tm.showPC(ax);
    catch e
        fprintf('%s: %s\n', newList{iExp}, e.message)
    end
end


end


%%
function list = get2pFiles(ExpRef)
end

function list = getTimelineFiles(ExpRef)
end

function list = getExpInfoFiles(ExpRef)
end

function list = getHardwareInfoFiles(ExpRef)
end

function list = getTrodesFiles(ExpRef)
end

function list = getEyeTrackingFiles(ExpRef)
end