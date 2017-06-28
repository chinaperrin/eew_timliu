function [pathName,recordName,restName] = splitFullName(fullName)

% pathName
slash_idx = regexp(fullName,'/');
traceName = fullName(slash_idx(end)+1:end);
pathName  = fullName(1:slash_idx(end));


[traceOrigin] = find_traceOrigin(fullName);

% JAPAN - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (strcmp(traceOrigin,'japan'))
    
    ptIdx      = regexp(traceName,'\.');
    recordName = traceName(1:ptIdx(1)-1);
    restName   = traceName(ptIdx(1):end);

    
% CALIFORNIA  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif (strcmp(traceOrigin,'cali'))
    
    recordName = traceName(1:end-5);
    restName   = traceName(end-4:end);
    
% NGA   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif (strcmp(traceOrigin,'nga'))
    
     uscIdx     = regexp(traceName,'_');
     recordName = traceName(1:uscIdx(end));
     restName   = traceName(uscIdx(end)+1:end);
end