function [neic_idx,differences] = find_neic_list_entry(eqName,someFileName,neicList)


% INPUT
% eqName        = name of eq-directory
% someFileName  = file name of some random waveform file of that eq (for reading header)
% neicList      = NEIC file list for reading correct magnitudes
% maxDist_time  = maximum distance in time for two list entries to be considered the same event
% maxDist_space = maximum distance in space for two list entries to be considered the same event


neicLat = neicList.eqLat;
neicLon = neicList.eqLon;
neicM   = neicList.m;

% Read some knt waveform file header from event 
[~,meta] = read_ascii_wform_jp(someFileName,17,1);
kntLat   = meta.eqLat;
kntLon   = meta.eqLon;
kntM     = meta.m;
kntTime  = meta.eqTime;


%% DIFFERENCES

% NEIC dates in seconds
n         = numel(neicList.eqLat);
neicDates = zeros(n,1);
for i = 1:n
    
    % Compute all neic time in serial dates / seconds  
    seconds      = neicList.time{i};
    seconds      = str2double(seconds(1:2))*3600 + str2double(seconds(4:5))*60 + str2double(seconds(7:9));
    neicDates(i) = datenum(neicList.date(i))*24*3600 + seconds; 
end

neicDates = neicDates + 9*3600;            % Time difference (as calibrated 
                                           % with the Tohoku event)

% k-net date in seconds
yr    = eqName(1:4);
mt    = eqName(5:6);
dy    = eqName(7:8);
hr    = eqName(9:10);
mn    = eqName(11:12);
sc    = eqName(13:14);
kDate = datenum([yr,'/',mt,'/',dy,',',hr,':',mn,':',sc])*24*3600;


dt   = kDate-neicDates; % [seconds]
dlat = kntLat-neicLat;
dlon = kntLon-neicLon;
dm   = kntM-neicM; 



% Identify events within +/- 50km (~ 0.5deg) & those which happened on the same hour/day
idx_near = find( (abs(dlat)<0.5)   & (abs(dlon)<0.5) );
idx_1h   = find( (abs(dt)<3600)    & (abs(dlat)<0.5) & (abs(dlon)<0.5) );
idx_1d   = find( (abs(dt)<3600*24) & (abs(dlat)<0.5) & (abs(dlon)<0.5) );



% FIGURE   -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
figure(889); clf; hold on
latmin=30; latmax=45; lonmin=128; lonmax=145;
hj = plot_japan_map(latmin,latmax,lonmin,lonmax);
plot(neicLon,neicLat,'.r','markerSize',5)
plot(kntLon,kntLat,'pk','markerFaceColor','y','markerSize',14)
title(['M ',num2str(kntM),' (knt-magnitude)'],'fontSize',15)

plot(neicLon(idx_near),neicLat(idx_near),'ob','markerSize',14)
axis([kntLon-0.7 kntLon+0.7 kntLat-0.7 kntLat+0.7])

if (~isempty(idx_1h))
    plot(neicLon(idx_1d),neicLat(idx_1d),'^k','markerSize',14,'lineWidth',2)
    plot(neicLon(idx_1h),neicLat(idx_1h),'vr','markerSize',14,'lineWidth',2)
    
    
    % Choose differences
    [val,idx] = min(abs(dt(idx_1h)));
    neic_idx  = idx_1h(idx);
    plot(neicLon(neic_idx),neicLat(neic_idx),'vk','markerFaceColor','y','markerSize',10)
    
    
    fprintf(1,['Found ',num2str(numel(idx_1h)),' events which happened within 50km and 1hr from M',num2str(kntM),' target event. Differences are:\n'])
    
    lat_offset=0;
    lon_offset=0.1;
    for i = 1:numel(idx_1h)
        strng = strcat(['dT = ',num2str(dt(idx_1h(i))),'  /  dM = ',num2str(dm(idx_1h(i)))]);
        text(neicLon(idx_1h(i))+lon_offset,neicLat(idx_1h(i))+lat_offset,strng,'fontSize',15)
        fprintf(1,['    dT = ',num2str(dt(idx_1h(i))),'  /  dM = ',num2str(dm(idx_1h(i))),'\n'])
    end
    
    
    differences.dt   = dt(idx_1h);
    differences.dm   = dm(idx_1h);
    differences.dlat = dlat(idx_1h);
    differences.dlon = dlon(idx_1h);
   
else
    
    fprintf(1,['Found no event which happened within 50km and 1hr from M',num2str(kntM),' target event\n'])

    neic_idx         = 99999;
    differences.dt   = 99999;
    differences.dm   = 99999;
    differences.dlat = 99999;
    differences.dlon = 99999;
end


%% Appendix
% dt = abs(kDate - neicDates);
% [dt,neic_idx] = min(dt_dates);
% fprintf(1,'How about doubles?\n')
%     
% if (dt>maxDist_t0)
%     neic_idx = [];
%     fprintf(1,['Closest event in neic-list is more than ',num2str(maxDist_t0),'sec apart\n'])
% else
%     1+1;
% end