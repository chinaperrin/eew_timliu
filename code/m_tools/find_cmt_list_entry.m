function [out] = find_cmt_list_entry(CMT,eqLat,eqLon,eqM,eqDate,eqTime,tPrime,mPrime,latPrime,lonPrime)

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
cmtDatString = strcat(CMT.date,',',CMT.t0);
cmtDatSec    = datenum(cmtDatString)*24*3600;

% Differences between target trace and all entries of cmt list
dt   = abs(eqDatSec - cmtDatSec);
dm   = abs(eqM      - CMT.Mw);
dlat = abs(eqLat    - CMT.lat);
dlon = abs(eqLon    - CMT.lon);

idx = find((dt<tPrime) & (dm<mPrime) & (dlat<latPrime) & (dlon<lonPrime));
 
% Write info into out structure
out.listIdx = idx;

if (numel(idx)==1)
    out.date  = CMT.date{idx};
    out.t0    = CMT.t0{idx};
    out.M0    = CMT.M0(idx);
    out.Mw    = CMT.Mw(idx);
    out.mb    = CMT.mb(idx);
    out.MS    = CMT.MS(idx);
    out.lat   = CMT.lat(idx);
    out.lon   = CMT.lon(idx);
    out.depth = CMT.depth(idx);
    
    if ~o.quiet
        %fprintf(1,'\n----------------------------------------------\n')
        %fprintf(1,['kNet-header:\t\t ',eqDate,' ',eqTime,'\n'])
        fprintf(1,['     found cmt entries:\t\t',out.date,' ',out.t0,' Mw_cnt ',num2str(out.Mw,'%3.1f'),'\n'])
    end

else
    fprintf(1,'     found cmt entries:\t\t --\n')
    %fprintf(1,['None or several corresponding events found in cmt list:\n'])
    fprintf(1,num2str(idx'))
end