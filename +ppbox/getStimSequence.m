function stimSequence=getStimSequence(animal, iseries, iexp)

SetDefaultDirs;

if nargin == 1
    if isstruct(animal)
        info = animal;
        animal=info.subject;
        iseries=str2num(info.expDate(info.expDate~='-'));
        iexp=info.exp;
    else
        expRef = animal;
        [animal, dd, iexp] = dat.parseExpRef(expRef);
        iseries = str2num(datestr(dd, 'yyyymmdd'));
    end    
elseif nargin == 3
    % do nothing
end
p=ProtocolLoad(animal, iseries, iexp);

stimSequence.labels=cell(p.nstim, 1);
stimSequence.seq=zeros(p.nrepeats*p.nstim, 1);

% % building the labels (according to the active parameters)
nActivePars=length(p.activepars);
for iStim=1:p.nstim
    for iActivePar=1:nActivePars
        ind=p.activepars{iActivePar};
        ind=(ind(min(length(ind), 2)));
        if iActivePar==1
            stimSequence.labels{iStim}=sprintf('%s = %d', p.parnames{ind}, p.pars(ind, iStim));
        else
            stimSequence.labels{iStim}=sprintf('%s, %s = %d', stimSequence.labels{iStim}, p.parnames{ind}, p.pars(ind, iStim));
        end
    end
    if ismember(iStim, p.blankstims)
        stimSequence.labels{iStim}='blank';
    end
end

% building the labels (for x1 and y1, as required by matteo's code)
% p.activepars={8, 9};
% nActivePars=length(p.activepars);
% for iStim=1:p.nstim
%     for iActivePar=1:nActivePars
%         ind=p.activepars{iActivePar};
%         if iActivePar==1
%             stimSequence.labels{iStim}=sprintf('%d) x1 = %d', iStim, p.pars(ind, iStim));
%         else
%             stimSequence.labels{iStim}=sprintf('%s\ny1 = %d', stimSequence.labels{iStim}, p.pars(ind, iStim));
%         end
%     end
% end

[ss ii]=sort(p.seqnums(:));
stimSequence.seq=repmat([1:p.nstim]', p.nrepeats, 1);
stimSequence.seq=stimSequence.seq(ii);


            
