function [out] = find_fnet_list_entry(FNET,eqLat,eqLon,eqM,eqDate,eqTime,tPrime,mPrime,latPrime,lonPrime)

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

% Eq info from fNet event list
fDatSec    = datenum(FNET.date)*24*3600;


% Differences between target trace and all entries of fNet list
dt   = abs(eqDatSec - fDatSec);
dm   = abs(eqM      - FNET.jmaM );
dlat = abs(eqLat    - FNET.eqLat);
dlon = abs(eqLon    - FNET.eqLon);

idx = find((dt<tPrime) & (dm<mPrime) & (dlat<latPrime) & (dlon<lonPrime));

% Write fNet info into out structure
out.listIdx = idx;

if (numel(idx)==1) 
    
    out.date   = FNET.date{idx};
    out.jmaM   = FNET.jmaM(idx);
    out.lat    = FNET.eqLat(idx);
    out.lon    = FNET.eqLon(idx);
    out.jmaZ   = FNET.jmaZ(idx);
    
    out.strike = FNET.strike(idx);
    out.dip    = FNET.dip(idx);
    out.rake   = FNET.rake(idx);

    out.Mw     = FNET.Mw(idx);
    out.M0     = FNET.M0(idx);
    out.MTZ    = FNET.MTZ(idx);
    
    if ~o.quiet
        %fprintf(1,'\n----------------------------------------------\n')
        %fprintf(1,['kNet-header:\t\t ',eqDate,' ',eqTime,'\n'])
        fprintf(1,['     found fnet entries:\t',strrep(out.date,',',' '),' Mw_fnt ',num2str(out.Mw),'\n'])
    end
else
    fprintf(1,'     found fnet entries:\t --\n')
    %fprintf(1,['None or several corresponding events found in fNet list:\n'])
    fprintf(1,num2str(idx'))
end