function ipick = get_miniList_indices(trList,ntr_per_dataSet)

% Make a mini list that includes all data sets
dsnames = unique(trList.dataSetName);
nds     = numel(dsnames);
ipick   = [];
for ids = 1:nds
    itmp = find(strcmp(trList.dataSetName,dsnames{ids}));
    if numel(itmp)>=ntr_per_dataSet; ipick = [ipick; itmp(1:ntr_per_dataSet)];
    else                             ipick = [ipick; itmp(1:numel(itmp))];
    end
end

% Add data from non-90ies SCSN
idxSCSN = find(strcmp(trList.dataSetName,'scsn'));
nscsn   = numel(idxSCSN);
irand   = randi(nscsn,ntr_per_dataSet,1);
itmp    = idxSCSN(irand);
ipick   = [ipick; itmp];