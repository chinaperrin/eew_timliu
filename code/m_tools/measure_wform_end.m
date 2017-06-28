%function [endIndex] = measure_wform_end(trFullName,trList)

%%
%  HOW TO DEFINE WAVEFORM END?
%  FOR LARGE EVENTS, SEISMOGRAM RINGS FOREVER
%%


o.debug = 1;

jpIdx  = find(strcmp(zList.dataSetName,'kNet'));
jpList = zList.selectSubList(jpIdx);
idx_nf = find(jpList.hypDist<30);
jpList = jpList.selectSubList(idx_nf);
% Find trList-index
itr = 111;
zFullName = jpList.fullName{itr};

% Load raw trace
[zS_raw,meta] = read_any_trace_noproc(zFullName);
zS_raw        = zS_raw-mean(zS_raw);

% Load processed trace
[zS,~,meta] = read_any_trace_proc(zFullName,jpList);
zPxIdx      = meta.ppxIdx;
tz          = meta.t;
sr          = meta.sr;
dt          = 1/sr;
ztpx        = tz(zPxIdx);

if o.debug
    figure(111); clf; hold on; grid on; box on;
    plot(tz,zS_raw,'b')
    plot(tz,zS.acc,':k','lineWidth',1)
    legend('raw','processed')
end

sCum   = cumsum(abs(zS_raw));
ds     = diff(sCum);                 % Change in cumulative absolute signal
ds_med = median(ds(1:zPxIdx));       % Median change before p-onset
tdiff  = tz(1:end-1)+dt;

figure(112); clf;
subplot(3,1,1); hold on; grid on; box on;
plot(tz,zS.vel)

subplot(3,1,2); hold on; grid on; box on;
plot(tz,sCum)

subplot(3,1,3); hold on; grid on; box on;
plot(tdiff,ds,'k')
line([tdiff(1) tdiff(zPxIdx)],     [ds_med ds_med],'color','r')
line([tdiff(end-1000) tdiff(end)], [10*ds_med 10*ds_med],'color','r')

% Measure noise

% Compare to what is saved in trList

% Make running window over wform, start at ts + (ts-tp)