function [out] = find_cmt_list_entry(cmtListFullName,eqLat,eqLon,eqM,eqDate,eqTime,tPrime,mPrime,latPrime,lonPrime)

o.quiet = false;

% Set similarity thresholds
if nargin==6
    tPrime   = 60;      % [sec]
    mPrime   = .5;      % [m-units]
    latPrime = .5;      % [deg]
    lonPrime = .5;      % [deg]
end

eqDatString = strcat(eqDate,',',eqTime);
eqDatSec    = datenum(eqDatString)*24*3600;
eqDatSec    = eqDatSec - 9*3600;                % Change from JST to UTC

% Eq info from cmt event list
[cmtList]    = import_cmtList(cmtListFullName);
cmtDatString = strcat(cmtList.date,',',cmtList.time);
cmtDatSec    = datenum(cmtDatString)*24*3600;

% Differences between target trace and all entries of cmt list
dt   = abs(eqDatSec - cmtDatSec);
dm   = abs(eqM      - cmtList.mb);
dlat = abs(eqLat    - cmtList.eqLat);
dlon = abs(eqLon    - cmtList.eqLon);

idx = find((dt<tPrime) & (dm<mPrime) & (dlat<latPrime) & (dlon<lonPrime));

% Write info into out structure
out.listIdx = idx;

if (numel(idx)==1)
    out.date    = cmtList.date{idx};
    out.t0      = cmtList.time{idx};
    out.mb      = cmtList.mb(idx);
    out.eqLat   = cmtList.eqLat(idx);
    out.eqLon   = cmtList.eqLon(idx);
    out.eqZ     = cmtList.eqZ(idx);
    
    if ~o.quiet
        fprintf(1,'\n----------------------------------------------\n')
        fprintf(1,['kNet-header:\t\t ',eqDate,' ',eqTime,'\n'])
        fprintf(1,['found cmt entries:\t',out.date,' ',out.t0,'\n'])
    end

else
    fprintf(1,['None or several corresponding events found in cmt list:\n'])
    fprintf(1,num2str(idx'))
end