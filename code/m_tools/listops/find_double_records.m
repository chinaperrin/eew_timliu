function doubleIdxList = find_double_records(TraceList)

%thisList = TraceList.selectSubList(1:1500);

% Single out record names (without path)
slashArray     = cellfun(@(x) regexp(x,'/'), TraceList.fullName,'uniformOutput',0);
lastSlashIdx   = cell2mat(cellfun(@(x) x(end), slashArray,'uniformOutput',0));
%allRecordNames = cellfun(@(x) x(lastSlashIdx+1:end), TraceList.fullName,'uniformOutput',0); didnt work right, dunno why

ntr            = numel(TraceList.eq.m);
allRecordNames = cell(ntr,1);
for itr = 1:ntr
    allRecordNames{itr} = TraceList.fullName{itr}(lastSlashIdx(itr)+1:end);
end


% Single out event names
eqNameList = unique(TraceList.eq.name);
neq        = numel(eqNameList);
doubleIdxList = [];
for ieq = 1:neq
    thisEqName = eqNameList{ieq};
    idx = find(strcmp(thisEqName,TraceList.eq.name));
    
    recordNames = allRecordNames(idx);
    nrecs       = numel(recordNames);
    alreadyProc = false(nrecs,1);
    irec        = 1;
    
    while irec<nrecs
        if ~alreadyProc(irec)
            thisName = recordNames(irec);
            idxMatch = find(strcmp(thisName,recordNames));
            if numel(idxMatch)>1
                doubleIdxList = [doubleIdxList; idx(idxMatch(2:end))];
                %check: TraceList.fullName(idx(idxMatch))
            end
            alreadyProc(idxMatch) = true;
        end
        irec=irec+1;
    end
end