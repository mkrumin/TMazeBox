function stimMatrix=buildStimMatrix(stimSequence, stimTimes, timeAxis)

nSamples=length(timeAxis);
nStim=length(stimSequence.labels);

stimMatrix=false(nStim, nSamples);

for iStim=1:length(stimSequence.seq)
    ind=(timeAxis>=stimTimes.onset(iStim) & timeAxis<=stimTimes.offset(iStim));
    stimMatrix(stimSequence.seq(iStim), ind)=true;
end



