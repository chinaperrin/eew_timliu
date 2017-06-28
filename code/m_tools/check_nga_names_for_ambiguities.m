function [c_ambiguous] = check_nga_names_for_ambiguities(TraceList)

fprintf(1,'Going through all TraceList entries and counting the number of corecords ...\n\n')
c_ambiguous = 0;

ntr = size(TraceList.m,1);
% For all traces in list ...
for itr = 1:ntr
    
    fprintf(1,['\n',num2str(itr),'/',num2str(ntr),' ...\t'])

    traceFullName = TraceList.fullName{itr};
    [recordName]  = get_recordName(traceFullName);
    
    % Identify corecorded traces
    hits     = regexp(TraceList.fullName,recordName);
    idxCorec = find(cellfun(@(x) ~isempty(x), hits));
    nCorec   = numel(idxCorec);
    if (nCorec~=3)
        c_ambiguous = c_ambiguous + 1;
        {TraceList.fullName{idxCorec}}'
        1+1;
    end
end