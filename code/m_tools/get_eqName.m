function [eqName] = get_eqName(traceFullName)


[traceOrigin] = find_traceOrigin(traceFullName);

% JAPAN - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (strcmp(traceOrigin,'japan'))
    
    slashIdx = regexp(traceFullName,'/');
    ptIdx    = regexp(traceFullName,'\.');
    eqName   = traceFullName(slashIdx(end)+7:ptIdx-1);
    
    
% CALIFORNIA  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif (strcmp(traceOrigin,'cali'))
    
    slashIdx = regexp(traceFullName,'/');
    ptIdx    = regexp(traceFullName,'\.');
    eqName   = traceFullName(slashIdx(end)+1:ptIdx(1)-1);
    
    
% NGA   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif (strcmp(traceOrigin,'nga'))
    
    slashIdx = regexp(traceFullName,'/');
    ptIdx    = regexp(traceFullName,'\.');
    eqName   = traceFullName(slashIdx(end)+1:ptIdx(1)-1);
    
end