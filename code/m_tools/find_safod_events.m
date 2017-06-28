clear all
addpath(genpath('~/programs/matlab/'))

DLAT = .03;     % Find events within +/- DLAT of borehole
DLON = .03;
DT   = 100;      % Search trigger list entries with start time +/-DT of catalog origin time

% Load triggerlist and catalog
listFullName     = '/scratch/memeier/data/nocal/safod/PH_triggers.txt';
triggerList      = import_SAFOD_PH_triggerList(listFullName);
[trStartDates,~] = jl2normaldate(triggerList.startDate*1000,'dd-mmm-yyyy');
trStartDates     = cellstr(trStartDates);
[trEndDates,~  ] = jl2normaldate(triggerList.endDate*1000,'dd-mmm-yyyy');
trEndDates       = cellstr(trEndDates);
%[normalDate status] = jl2normaldate(2009.256)
%[normalDate status] = jl2normaldate(triggerList.startDate(4404)*1000);


[ctlg] = import_NCAeqDD('/scratch/memeier/data/nocal/waldhauser/NCAeqDD.v201112.1.txt');





%% TIME WINDOW
% Time window of triggers in datenum format
dynum_start = datenum(trStartDates(1));
dynum_end   = datenum(trStartDates(end));


ctlg_dynum = datenum(ctlg.yr,ctlg.mt,ctlg.dy,ctlg.hr,ctlg.min,ctlg.sec);
lgc_time   = (ctlg_dynum>=dynum_start & ctlg_dynum<=dynum_end);


%% SPATIAL WINDOW
%lat_mh =   35.9742;
%lon_mh = -120.5521;
lat_ph =   35.9743;
lon_ph = -120.5521;

latmax = lat_ph + DLAT;
latmin = lat_ph - DLAT;
lonmax = lon_ph + DLON;
lonmin = lon_ph - DLON;

lgc_space = (ctlg.lat>=latmin & ctlg.lat<=latmax & ctlg.lon>=lonmin & ctlg.lon<=lonmax);


%% Extract sub-catalog
idx     = find(lgc_time & lgc_space);
sfdCtlg = ctlg(idx,:);
numel(idx);

[val,idx] = sort(sfdCtlg.m,'descend');
sfdCtlg   = sfdCtlg(idx,:);

idx2 = find(sfdCtlg.m>=2); 
idx3 = find(sfdCtlg.m>=3 & sfdCtlg.m<4); 
idx4 = find(sfdCtlg.m>=4 & sfdCtlg.m<5);


%% Plot map
figure(333); hold on;
plot_Cali_map
plot(sfdCtlg.lon,sfdCtlg.lat,'xk')
plot(sfdCtlg.lon(idx2),sfdCtlg.lat(idx2),'pk','markerSize',14,'markerFaceColor','m')
plot(sfdCtlg.lon(idx3),sfdCtlg.lat(idx3),'sk','markerSize',16,'markerFaceColor','m')
plot(sfdCtlg.lon(idx4),sfdCtlg.lat(idx4),'vk','markerSize',18,'markerFaceColor','m')
plot(lon_ph,lat_ph,'ok','markerFaceColor','r','markerSize',10)
%plot(lon_mh,lat_mh,'dk','markerFaceColor','b','markerSize',10)
set(gca,'xlim',[lonmin lonmax],'ylim',[latmin latmax])




% Largest events sufficiently close to borehole
m2ctlg = sfdCtlg(idx2,:);
neq    = size(m2ctlg,1);

% Find them in triggerList
tstart       = cellfun(@(x,y) [x,' ',y],trStartDates,triggerList.startTime,'uniformOutput',0);
tstart_dynum = datenum(tstart);
tend         = cellfun(@(x,y) [x,' ',y],trEndDates,triggerList.endTime,'uniformOutput',0);
tend_dynum   = datenum(tend);

% Find Julian day without year
julian_day = juliandate([m2ctlg.yr,m2ctlg.mt,m2ctlg.dy])-juliandate([m2ctlg.yr,ones(neq,1),ones(neq,1)])+1;

for ieq = 1:neq
    ctlg_dynum = datenum(m2ctlg.yr(ieq),m2ctlg.mt(ieq),m2ctlg.dy(ieq),m2ctlg.hr(ieq),m2ctlg.min(ieq),m2ctlg.sec(ieq));
    dt         = ctlg_dynum-tstart_dynum;
    [val,idx]  = find(abs(dt)<DT); 
    whos idx
end


% No triggers in list for 2002.264
% ...
% 1176 67001176    2003.293,04:52:18.2165  2003.293,04:52:33.2160
% 1177 67001177    2003.293,04:52:33.2255  2003.293,04:52:48.2250
% 1178 67001178    2003.293,08:03:47.7690  2003.293,08:04:02.7685
% 1179 67001179    2003.293,09:09:37.3570  2003.293,09:09:52.3565
% 1180 67001180    2003.293,09:09:52.3715  2003.293,09:10:07.3710
% 1181 67001181    2003.293,10:31:40.1505  2003.293,10:31:55.1500
% 1182 67001182    2003.293,11:25:44.6145  2003.293,11:25:59.6140 --> dloaded
% 1183 67001183    2003.293,11:25:59.6290  2003.293,11:26:14.6285 --> dloaded
% 1184 67001184    2003.293,11:31:27.3280  2003.293,11:31:42.3275 --> dloaded
% 1185 67001185    2003.293,12:01:29.1525  2003.293,12:01:44.1520
% ...
% 1186 67001186    2003.294,00:19:14.0775  2003.294,00:19:29.0770
% 1187 67001187    2003.294,05:16:52.0085  2003.294,05:17:07.0080
% 1188 67001188    2003.294,08:30:43.8970  2003.294,08:30:58.8965
% 1189 67001189    2003.294,09:00:14.2570  2003.294,09:00:29.2565 --> dloaded
% 1190 67001190    2003.294,13:43:50.9510  2003.294,13:44:05.9505
% 1191 67001191    2003.294,15:53:16.9320  2003.294,15:53:31.9315
% 1192 67001192    2003.294,15:57:20.8145  2003.294,15:57:35.8140
% 1193 67001193    2003.294,23:42:58.7490  2003.294,23:43:13.7485




