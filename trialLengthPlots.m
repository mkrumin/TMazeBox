
folder = 'G:\DATA\';

aavDatasets = {...
    '2014-08-02_2023_MK012_TMwExtras.mat'; ...
    '2014-08-04_1827_MK012_TMwExtras.mat'; ...
    '2014-08-05_2228_MK012_TMwExtras.mat'; ...
    '2014-08-08_2251_MK012_TMwExtras.mat'; ...
    '2014-08-11_2133_MK012_TMwExtras.mat'; ...
    '2014-08-13_2023_MK012_TMwExtras.mat'; ...
    '2014-08-15_1931_MK012_TMwExtras.mat'; ...
    '2014-08-02_2203_MK014_TMwExtras.mat'; ...
    '2014-08-05_1937_MK014_TMwExtras.mat'; ...
    };

vglutDatasets = {...
    '2017-06-12_1420_JL005_TMwExtras.mat'; ...
    '2017-07-15_1708_JL008_TMwExtras.mat'; ...
    '2017-07-27_1433_JL008_TMwExtras.mat'; ...
    '2017-08-12_1056_JL008_TMwExtras.mat'; ...
    '2017-09-18_1707_JL008_TMwExtras.mat'; ...
    '2017-09-23_1539_JL008_TMwExtras.mat'; ...
    };

tripleDatasets = {...
    '2015-07-03_2127_MK020_TMwExtras.mat'; ...
    '2015-07-07_2127_MK020_TMwExtras.mat'; ...
    '2015-07-30_2010_MK020_TMwExtras.mat'; ...
    '2015-07-31_1851_MK020_TMwExtras.mat'; ...
    '2015-08-02_1709_MK020_TMwExtras.mat'; ...
    '2015-08-05_1930_MK020_TMwExtras.mat'; ...
    '2015-08-07_1821_MK020_TMwExtras.mat'; ...
    '2015-07-03_1632_MK022_TMwExtras.mat'; ...
    '2015-08-02_2051_MK022_TMwExtras.mat'; ...
    '2015-06-19_2155_MK023_TMwExtras.mat'; ...
    '2015-07-18_154_MK023_TMwExtras.mat';
    };

allPPC = [aavDatasets; tripleDatasets; vglutDatasets];

allFiles = dir(sprintf('%s*_TMwExtras.mat', folder));
files{1} = allFiles(ismember({allFiles.name}, aavDatasets));
files{2} = allFiles(ismember({allFiles.name}, tripleDatasets));
files{3} = allFiles(ismember({allFiles.name}, vglutDatasets));

% idx = ismember({allFiles.name}, allPPC);
% % idx = ismember({allFiles.name}, '2017-09-23_1539_JL008_TMwExtras.mat');
% % idx = ismember({allFiles.name}, '2015-08-02_1709_MK020_TMwExtras.mat');
% idx = ismember({allFiles.name}, '2014-08-05_1937_MK014_TMwExtras.mat');
%
% files = mat2cell(allFiles(idx), ones(size(allFiles(idx))), 1);

nGroups = length(files);

groupName = {'AAV', '3tg', 'VGlut', 'All mice'};
% for iGroup = 1:nGroups
%     groupName{iGroup} = files{iGroup}(1).name(1:21);
% end
groupName{nGroups+1} = 'All mice';
groups2plot = 1:length(groupName);

%% analysis

tic
res = cell(nGroups, 1);
for iGroup = 1:nGroups
    fprintf('group %1.0f/%1.0f\n', iGroup, nGroups);
    nFiles = length(files{iGroup});
    for iFile = 1:nFiles
        fprintf('file %1.0f/%1.0f\n', iFile, nFiles);
        try
            load(fullfile(folder, files{iGroup}(iFile).name));
        catch
            warning(sprintf('failed to load %s \n', files{iGroup}(iFile).name))
            continue
        end

        res{iGroup}(iFile) = getTrialDur(TM);

    end
    res{iGroup} = res{iGroup}(:);
end


toc

%% define some parameters
resFull = res;
resFull{end+1} = cell2mat(res);

allExpRefs = {resFull{end}.ExpRef}';
animalNames = {'MK012', 'MK014', 'MK020', 'MK022', 'MK023', 'JL005', 'JL008'};
genotype = {'AAV','AAV', '3g', '3g', '3g', 'vGlut','vGlut'};

animalIdx = cellfun(@contains, ...
    repmat(allExpRefs, 1, length(animalNames)),...
    repmat(animalNames, length(allExpRefs), 1));
% add column of ones as last 'animal' - for all animals
animalIdx = cat(2, animalIdx, ones(size(animalIdx, 1), 1));
animalNames = [animalNames, {'All'}];
genotype = [genotype, {'Mix'}];
nAnimals = length(animalNames);

%%
figure
nRows = 2;
nColumns = ceil(nAnimals/nRows);
for iAnimal = 1:nAnimals
    dur = cell2mat({resFull{end}(find(animalIdx(:,iAnimal))).dur}');
    isFinished = cell2mat({resFull{end}(find(animalIdx(:,iAnimal))).isFinished}');
    subplot(nRows, nColumns, iAnimal);
    histogram(dur(isFinished));
    medDur = median(dur(isFinished));
    madDur = mad(dur(isFinished), 1);
    title({animalNames{iAnimal}; sprintf('%3.1f\\pm%3.1f', medDur, madDur)})
end


