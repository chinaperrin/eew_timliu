function availList = get_available_traceList(trList)
% Go through trList and check if waveform files exist. If not, try to
% scp them to local machine. If that doesn't work, throw them out of list

ntr        = numel(trList.m);
throwMeOut = false(ntr,1);
for itr = 1:ntr
    
    print_iteration_numbers(itr,ntr,'tens')
    traceFullName = trList.fullName{itr};
    if ~exist(traceFullName,'file')
        
        status = scp_wform(traceFullName);
        if status~=0;
            fprintf(1,'not found: %s\n',traceFullName)
            throwMeOut(itr) = true;
        end
    end
end
idxKeep   = find(throwMeOut==0);
availList = trList.selectSubList(idxKeep);
