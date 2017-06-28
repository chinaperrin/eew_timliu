function [indices] = find_all_event_traces(traceFullName,trList)


[traceOrigin] = find_traceOrigin(traceFullName);

% JAPAN - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (strcmp(traceOrigin,'japan'))
    
    slashIdx = regexp(traceFullName,'/');
    ptIdx    = regexp(traceFullName,'\.');
    eqName   = traceFullName(slashIdx(end)+7:ptIdx-1);
    
    sameEvent = regexp(trList.fullName,eqName,'match');
    indices   = find(cellfun(@(x) ~isempty(x), sameEvent));
    
    
% CALIFORNIA  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif (strcmp(traceOrigin,'cali'))
    
    slashIdx = regexp(traceFullName,'/');
    ptIdx    = regexp(traceFullName,'\.');
    eqName   = traceFullName(slashIdx(end)+1:ptIdx(1)-1);
    
    sameEvent = regexp(trList.fullName,eqName,'match');
    indices   = find(cellfun(@(x) ~isempty(x), sameEvent));

% CALIFORNIA  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif (strcmp(traceOrigin,'nocal'))
    
    slashIdx = regexp(traceFullName,'/');
    ptIdx    = regexp(traceFullName,'\.');
    eqName   = traceFullName(slashIdx(end-1)+1:ptIdx(1)-7);
    
    sameEvent = regexp(trList.fullName,eqName,'match');
    indices   = find(cellfun(@(x) ~isempty(x), sameEvent));

    
% NGA   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif (strcmp(traceOrigin,'nga'))
    
    slashIdx = regexp(traceFullName,'/');
    uscIdx   = regexp(traceFullName,'\_');
    eqName   = traceFullName(slashIdx(end)+1:uscIdx(1)-1);
    
    sameEvent = regexp(trList.fullName,eqName,'match');
    indices   = find(cellfun(@(x) ~isempty(x), sameEvent));

% WENCHUAN    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif (strcmp(traceOrigin,'wenchuan'))
   
    % 20080512142804_Ms8p0_051WCW_ORIG_UD/NS/EW.dat
    slashIdx = regexp(traceFullName,'/');
    uscIdx   = regexp(traceFullName,'\_');
    eqName   = traceFullName(slashIdx(end)+1:uscIdx(2)-1);
    
    sameEvent = regexp(trList.fullName,eqName,'match');
    indices   = find(cellfun(@(x) ~isempty(x), sameEvent));
end