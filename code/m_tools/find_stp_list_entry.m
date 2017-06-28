function [stp_idx] = find_stp_list_entry(eqName,stpList,maxDist_t0)


% F-net dates in seconds
stpDates = datenum(stpList.date)*24*3600;
stpDates = stpDates + 9*3600;               % Time difference (as calibrated 
                                            % by the Tohoku event)

yr = eqName(1:4);
mt = eqName(5:6);
dy = eqName(7:8);
hr = eqName(9:10);
mn = eqName(11:12);
sc = eqName(13:14);

% k-net date in seconds
kDate = datenum([yr,'/',mt,'/',dy,',',hr,':',mn,':',sc])*24*3600;

dt_dates     = abs(kDate - stpDates);
[dt,stp_idx] = min(dt_dates);

    fprintf(1,'How about doubles?\n')
    
if (dt>maxDist_t0)
    stp_idx = [];
    fprintf(1,'Closest event in stp-list is more than 300sec apart\n')
else
    1+1
end

