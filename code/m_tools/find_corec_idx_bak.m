function [idx_Z,idx_E,idx_N] = find_corec_idx(traceFullName,traceList,includeOtherInstr)

% This function finds the traceList indices of the three (?) waveform files
% that were recorded at the same site with the same instrument (HH_ pattern). 
% If <includeOtherInstr> is set to true, also other instruments from the same 
% Site are included (H__ pattern).

if nargin==2; includeOtherInstr=0; end

idx_Z   = [];
idx_N   = [];
idx_E   = [];

% Find trace in traceList            -  -  -  -  -  -  -  -  -  -  -  -  -
%glIdx = find(strcmp(traceFullName,traceList.fullName));

% Find out what sort of waveform it is ------------------------------------
ptIdx         = regexp(traceFullName,'\.');
extension     = traceFullName(ptIdx(end):end);
jp_extensions = {'.EW','.EW1','.EW2','.NS','.NS1','.NS2','.UD','.UD1','.UD2'};


if strcmp(extension,'.sac')       % - - - - - - - - - - SOUTHERN CALIFORNIA

    % Find co-recorde traces
    slashIdx   = regexp(traceFullName,'/');
    if ~includeOtherInstr; recordName = traceFullName(slashIdx(end)+1:end-5);
    else                   recordName = traceFullName(slashIdx(end)+1:end-6);
    end
    isCorec    = regexp(traceList.fullName,recordName,'match');
    idxCorec   = find(cellfun(@(x) ~isempty(x), isCorec));
    
    instrument = recordName(end);
    
    % Find out which trace is which       
    names   = traceList.fullName(idxCorec);
    
    matches = regexp(names,'\Z.sac', 'end');                  
    idx_tmp = find(cellfun(@(x) ~isempty(x),matches));
    idx_Z   = idxCorec(idx_tmp);
    
    matches = regexp(names,'\N.sac', 'end');
    idx_tmp = find(cellfun(@(x) ~isempty(x),matches));
    idx_N   = idxCorec(idx_tmp);
    
    matches = regexp(names,'\E.sac', 'end'); 
    idx_tmp = find(cellfun(@(x) ~isempty(x),matches));
    idx_E   = idxCorec(idx_tmp);
    
    
elseif ismember(extension,jp_extensions)     % - - - - - - - - - - -  JAPAN
    % Find co-recorded traces
    %fprintf(1,'Finding co-recorded traces ... ')
    tmp        = regexp(traceFullName,'/','split');
    tmp2       = tmp{end};
    tmp3       = regexp(tmp2,'\.','split');
    recordName = tmp3{1};
    isCorec    = regexp(traceList.fullName,recordName);
    idxCorec   = find(cellfun(@(x) ~isempty(x), isCorec));
    
    % Find out which trace is which 
    names   = traceList.fullName(idxCorec);
    
    matches = regexp(names,'\.U', 'end');
    idx_tmp = find(cellfun(@(x) ~isempty(x),matches));
    idx_Z   = idxCorec(idx_tmp);
    
    matches = regexp(names,'\.N', 'end');
    idx_tmp = find(cellfun(@(x) ~isempty(x),matches));
    idx_N   = idxCorec(idx_tmp);
    
    matches = regexp(names,'\.E', 'end');
    idx_tmp = find(cellfun(@(x) ~isempty(x),matches));
    idx_E   = idxCorec(idx_tmp);

elseif strcmp(extension,'.AT2')             % - - - - - - - - - - NGA West1
    
    % Find co-recorded traces
    slash_idx     = regexp(traceFullName,'/');
    traceName     = traceFullName(slash_idx(end)+1:end);
    pt_idx        = regexp(traceName,'\.');
    recordName    = traceName(1:pt_idx(1));
    recordName    = strrep(recordName,'.','\.');
    
    isCorec    = regexp(traceList.fullName,recordName);
    idxCorec   = find(cellfun(@(x) ~isempty(x), isCorec));
    
    % Find out which trace is which 
    names   = traceList.fullName(idxCorec);
    
    matches = regexp(names,'\.V\.', 'end');
    idx_tmp = find(cellfun(@(x) ~isempty(x),matches));
   
    if (~isempty(idx_tmp))
        
        idx_Z   = idxCorec(idx_tmp);
        idxH    = idxCorec(idxCorec~=idx_Z);
        if (numel(idxH)==2)
            idx_N   = idxH(1);
            idx_E   = idxH(2);
        
        elseif (numel(idxH)==1)
            idx_N   = idxH(1);
            flg_noE = true;
        
        elseif (numel(idxH)==0)
            flg_noN = true;
            flg_noE = true;
        end
    end
    
else
    fprintf(1,'EXTENSION NOT FOUND, WHAT NOW?\n')
    pause
end    % ------------------------------------------------------------------    
