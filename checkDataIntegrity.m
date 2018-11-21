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

paperData = load('G:\Krumin_etal_2018_eLife.mat');
rootFolder = 'G:\Processing\';
nExps = length(allExpRefs);
for iExp = 1:nExps
    [subject, d, expNum] = dat.parseExpRef(allExpRefs{iExp});
    d = datestr(d, 'yyyy-mm-dd');
    folderName = fullfile(rootFolder, subject, d, num2str(expNum));
    files = dir(fullfile(folderName, '*_ROI.mat'));
    idx = cellfun(@isempty, (strfind({files.name}, 'rect')));
    files = files(idx);
    nCellsROI = NaN;
    nCellsS2P = NaN;
    nCellsPaper = NaN;
    for iFile = 1:length(files)
        load(fullfile(folderName, files(iFile).name));
        [~, idx] = ismember(meta.ROI.CellClasses, 's');
        nCellsROI = sum(~all(isnan(meta.F(:, find(idx)))));
        if isfield(meta, 'folderDat')
            try
                s2pdat = load(fullfile(meta.folderDat, meta.filenameDat));
            catch
                folderDat = strrep(meta.folderDat, '\F', '');
                s2pdat = load(fullfile(folderDat, meta.filenameDat));
            end
            nCellsS2P = sum([s2pdat.dat.stat.iscell]) - ...
                sum(all(isnan(s2pdat.dat.Fcell{1}(find([s2pdat.dat.stat.iscell]), :)), 2));
        end
        nCellsPaper = [paperData.data(iExp).traces(iFile).nCells];
        fprintf('%s:\t %g\t%g\t[%s]\n', allExpRefs{iExp}, nCellsROI, nCellsS2P, num2str(nCellsPaper));
    end
end