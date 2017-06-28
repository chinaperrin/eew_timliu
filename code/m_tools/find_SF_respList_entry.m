function [idx_resp] = find_SF_respList_entry(traceFullName,respList_SF)


%% Parse wform file name
ptIdx  = regexp(traceFullName,'\.');
eqYr   = str2double(traceFullName(ptIdx(1) -4:ptIdx(1 )-1));
eqJday = str2double(traceFullName(ptIdx(1) +1:ptIdx(2 )-1));
eqHr   = str2double(traceFullName(ptIdx(2) +1:ptIdx(3 )-1));
eqMn   = str2double(traceFullName(ptIdx(3) +1:ptIdx(4 )-1));
eqSec  = str2double(traceFullName(ptIdx(4) +1:ptIdx(6 )-1));
net    = traceFullName(ptIdx(6) +1:ptIdx(7 )-1);
sta    = traceFullName(ptIdx(7) +1:ptIdx(8 )-1);
loc    = traceFullName(ptIdx(8) +1:ptIdx(9 )-1);
cha    = traceFullName(ptIdx(9) +1:ptIdx(10)-1);
dunno  = traceFullName(ptIdx(10)+1:ptIdx(11)-1);
ext    = traceFullName(ptIdx(11)+1:end);

% Serial Date Number of record date as constructed from wform file name
% Test record: 2003.293.11.25.44.6145.SF.PH001.01.GP1.D.SA
% jdate = 2003293.4762108161
jdate       = eqYr*1e3+eqJday + 1/86400*(eqHr*3600+eqMn*60+eqSec);
[recDate,~] = jl2normaldate(jdate);
recDate_sdn = datenum(recDate);
fprintf(1,'8UNG: time computed from file name does not exactly match time in sac-header\n')
% seems to be somewhat off from the time stored in sac-header,
% on the order of single seconds --> enough for choosing right
% respList entry, but might give rise to problems in other situations




%% Find sub-list for appropriate channel
idxList = find(strcmp(net,respList_SF.net) &strcmp(sta,respList_SF.sta) &strcmp(cha,respList_SF.cha) &strcmp(loc,respList_SF.loc));
nresp   = numel(idxList);
%respList_SF(end-5:end,:)
%respList_SF(idx,:)
%respList_SF.ondate


%% In sub-list, find entry of right time period
% Serial date numbers of onset dates in respList entries
% Test: choose as ondate the same as recDate, but in format of date given in respList:
%       ondate= '2003.293.11:25:43'
onDates_sdn  = zeros(nresp,1);
offDates_sdn = zeros(nresp,1);
for iresp = 1:nresp
    
    ondate              = respList_SF.ondate{idxList(iresp)};
    j_ondate            = str2double(ondate(1:4))*1e3+str2double(ondate(6:8)) + 1/86400*( str2double(ondate(10:11))*3600 + str2double(ondate(13:14))*60 + str2double(ondate(16:17)));
    [ondate_string,~]   = jl2normaldate(j_ondate);
    onDates_sdn(iresp)  = datenum(ondate_string);
    %Test: onDate_sdn-recDate_sdn
    
    offdate             = respList_SF.offdate{idxList(iresp)};
    j_offdate           = str2double(offdate(1:4))*1e3+str2double(offdate(6:8)) + 1/86400*( str2double(offdate(10:11))*3600 + str2double(offdate(13:14))*60 + str2double(offdate(16:17)));
    [offdate_string,~]  = jl2normaldate(j_offdate);
    offDates_sdn(iresp) = datenum(offdate_string);
end

idx = find(recDate_sdn>onDates_sdn & recDate_sdn<offDates_sdn);
if ~isempty(idx); idx_resp = idxList(idx);
else              idx_resp = [];
    fprintf(1,'no matching respList entry found.\n')
end