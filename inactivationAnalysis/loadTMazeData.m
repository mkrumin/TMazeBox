function dataOut = loadTMazeData(files)

dataOut = struct('EXP', [], 'SESSION', []);
iData = 0;
for iExp = 1:length(files)
    data = load(files{iExp});
    % only take experiments with inactivations
    if data.EXP.optiStim
        iData = iData + 1;
        dataOut(iData) = data;
    end
end

dataOut = dataOut(:);
