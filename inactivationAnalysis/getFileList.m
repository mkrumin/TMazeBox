function fnList = getFileList(ls)

fnList = {};
for iList = 1:length(ls)

    [ExpRef, expDate] = dat.listExps(ls(iList).animalName);
    fn = dat.expFilePath(ExpRef, 'tmaze', 'master');
    
    [expDate, idx] = sort(expDate, 'ascend');
    ExpRef = ExpRef(idx);
    fn = fn(idx);
    
    % optogenetic exps didn't start before that
    idxValid = expDate >= datenum(ls(iList).startDate) & expDate <= datenum(ls(iList).endDate);
    if ~isempty(ls(iList).excludeDate)
        % exclude certain dates dates
        idxValid = idxValid & ~ismember(expDate, datenum(ls(iList).excludeDate));
    end
    
    if ~isempty(ls(iList).excludeSession)
        % exclude certain sessions
        idxValid = idxValid & ~ismember(ExpRef, ls(iList).excludeSession);
    end
    
%     idxValid = find(idxValid);
    
    fnList = cat(1, fnList, fn(idxValid));
    
end
