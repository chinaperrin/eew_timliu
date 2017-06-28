function [fnet_idx,differences] = find_fnet_entry(eqName,someFileName,fnetList,o_plot)


fnetLat = fnetList.eqLat;
fnetLon = fnetList.eqLon;
fnetM   = fnetList.jmaM;

% Read some knt waveform file header from event 
[~,meta] = read_ascii_wform_jp(someFileName,17,1);
kntLat   = meta.eqLat;
kntLon   = meta.eqLon;
kntM     = meta.m;
kntTime  = meta.eqTime;


% F-net dates in seconds
fDates = datenum(fnetList.date)*24*3600;

yr = eqName(1:4);
mt = eqName(5:6);
dy = eqName(7:8);
hr = eqName(9:10);
mn = eqName(11:12);
sc = eqName(13:14);

% k-net date in seconds
kDate = datenum([yr,'/',mt,'/',dy,',',hr,':',mn,':',sc])*24*3600;


dt   = kDate-fDates; % [seconds]
dlat = kntLat-fnetLat;
dlon = kntLon-fnetLon;
dm   = kntM-fnetM; 



% Identify events within +/- 50km (~ 0.5deg) & those which happened on the same hour/day
idx_near = find( (abs(dlat)<0.5)   & (abs(dlon)<0.5) );
idx_1h   = find( (abs(dt)<3600)    & (abs(dlat)<0.5) & (abs(dlon)<0.5) );
idx_1d   = find( (abs(dt)<3600*24) & (abs(dlat)<0.5) & (abs(dlon)<0.5) );


% FIGURE   -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
if (o_plot)
    figure(889); clf; hold on
    latmin=30; latmax=45; lonmin=128; lonmax=145;
    plot_japan_map(latmin,latmax,lonmin,lonmax);
    plot(fnetLon,fnetLat,'.r','markerSize',5)
    plot(kntLon,kntLat,'pk','markerFaceColor','y','markerSize',14)
    title(['Finding fNet-entry: M ',num2str(kntM),' (knt-magnitude)'],'fontSize',15)
    
    plot(fnetLon(idx_near),fnetLat(idx_near),'ob','markerSize',14)
    axis([kntLon-0.7 kntLon+0.7 kntLat-0.7 kntLat+0.7])
end

if (~isempty(idx_1h))
    
    % Choose differences
    [val,idx] = min(abs(dt(idx_1h)));
    fnet_idx  = idx_1h(idx);

    if(o_plot)
        plot(fnetLon(idx_1d),fnetLat(idx_1d),'^k','markerSize',14,'lineWidth',2)
        plot(fnetLon(idx_1h),fnetLat(idx_1h),'vr','markerSize',14,'lineWidth',2)
    
        plot(fnetLon(fnet_idx),fnetLat(fnet_idx),'vk','markerFaceColor','y','markerSize',10)
    end
    
    fprintf(1,['Found ',num2str(numel(idx_1h)),' events which happened within 50km and 1hr from M',num2str(kntM),' target event. Differences are:\n'])
    
    lat_offset=0;
    lon_offset=0.1;
    for i = 1:numel(idx_1h)
        strng = strcat(['dT = ',num2str(dt(idx_1h(i))),'  /  dM = ',num2str(dm(idx_1h(i)))]);
        text(fnetLon(idx_1h(i))+lon_offset,fnetLat(idx_1h(i))+lat_offset,strng,'fontSize',15)
        fprintf(1,['    dT = ',num2str(dt(idx_1h(i))),'  /  dM = ',num2str(dm(idx_1h(i))),'\n'])
    end
    
    
    differences.dt   = dt(idx_1h);
    differences.dm   = dm(idx_1h);
    differences.dlat = dlat(idx_1h);
    differences.dlon = dlon(idx_1h);
   
else
    fprintf(1,['Found no event which happened within 50km and 1hr from M',num2str(kntM),' target event\n'])
    
    fnet_idx         = [];
    differences.dt   = 99999;
    differences.dm   = 99999;
    differences.dlat = 99999;
    differences.dlon = 99999;
end