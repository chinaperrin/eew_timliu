function [blackIdxList,whiteIdxList,blComments] = find_blackListed_traces(TraceList,blackList)

nbl        = numel(blackList.eq.m);
ntr        = numel(TraceList.eq.m);
isblack    = false(ntr,1);
blComments = {};

for ibl = 1:nbl
    print_iteration_numbers(ibl,nbl,'hundreds')
    blFullName    = blackList.fullName{ibl};
    m             = blackList.eq.m(ibl);
    idxCandidates = find(TraceList.eq.m==m);
    
    if ~isempty(idxCandidates); idxBlackListed = find(cellfun(@(x) ~isempty(x), regexp(blFullName,TraceList.fullName(idxCandidates))));
    else                        idxBlackListed = [];
    end
    
    if ~isempty(idxBlackListed)
    
        % Identify corecords of black-listed traces
        idxInTraceList     = idxCandidates(idxBlackListed);         % idxBlackListed is relative to shorter list that is selected using idxCandidates
        traceFullName      = TraceList.fullName{idxInTraceList};
        if ~strcmp(blFullName,traceFullName); fprintf(1,'Black-list entry and found entry dont match. check.\n'); pause; end
        recordName         = get_recordName(traceFullName);    
        hits               = regexp(TraceList.fullName,recordName);
        idxCorecs          = find(cellfun(@(x) ~isempty(x), hits));
        isblack(idxCorecs) = true;
        blComments         = [blComments; blackList.comment{ibl}; blackList.comment{ibl}; blackList.comment{ibl}];
        if numel(idxCorecs)~=3; 
            fprintf(1,'8ung: not three corecs found. check.\n'); 
            1+1;
        end
    end
end

blackIdxList = find( isblack);
whiteIdxList = find(~isblack);