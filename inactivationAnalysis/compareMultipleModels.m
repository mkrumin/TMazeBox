function out = compareMultipleModels(contrast, behavior, idx, models, options)

nGroups = length(models);

for iGroup = 1:nGroups - 1
    for jGroup = iGroup + 1:nGroups
        out(iGroup, jGroup) = ...
            compareModels(contrast, behavior, idx([iGroup, jGroup]), ...
            models([iGroup, jGroup]), options);
    end
end
