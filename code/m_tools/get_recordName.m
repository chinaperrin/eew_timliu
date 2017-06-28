function [recordName] = get_recordName(traceFullName)


[traceOrigin] = find_traceOrigin(traceFullName);

% JAPAN - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (strcmp(traceOrigin,'japan'))
    
    slashIdx   = regexp(traceFullName,'/');
    ptIdx      = regexp(traceFullName,'\.');
    recordName = traceFullName(slashIdx(end)+1:ptIdx(1)-1);
    
% CALIFORNIA  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif (strcmp(traceOrigin,'cali'))
    
    slashIdx   = regexp(traceFullName,'/');
    recordName = traceFullName(slashIdx(end)+1:end-5);
    
% NGA   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif (strcmp(traceOrigin,'nga'))
    
     slashIdx   = regexp(traceFullName,'/');
     uscIdx     = regexp(traceFullName,'_');
     recordName = traceFullName(slashIdx(end)+1:uscIdx(end));
     %recordName = traceFullName(slashIdx(end)+1:uscIdx(2));
end