function [hf,ax] = prepare_singleRecTdpa_figure(tdpaValues,tRange)

global ftSize

% Set up figure
hf = figure(12); clf; box on; grid on; whitebg('w')
ax.a1 = subplot(4,2,[1,3,5,7]);
ax.a2 = subplot(4,2,2);
ax.a3 = subplot(4,2,4);
ax.a4 = subplot(4,2,6);
ax.a5 = subplot(4,2,8);

pos1 = get(ax.a1,'position'); pos2 = get(ax.a2,'position'); pos3 = get(ax.a3,'position');
pos4 = get(ax.a4,'position'); pos5 = get(ax.a5,'position');

yB = pos1(2);
yT = pos1(2)+pos1(4);
dy = (yT-yB)/19;  % 4dy for each of 4 subplots, 2dy for one gap

pos5(4)=4*dy; pos4(4)=4*dy; pos3(4)=4*dy; pos2(4)=4*dy;
pos4(2)=yB+7*dy;
pos3(2)=yB+11*dy;
pos2(2)=yB+15*dy;

set(ax.a5,'position',pos5); set(ax.a4,'position',pos4)
set(ax.a3,'position',pos3); set(ax.a2,'position',pos2)


% Plot tdpa-curves and S-arrivals
colours = tdpaValues{1}.colours;
nr = numel(tdpaValues);
for ir = 1:nr
    t          = tdpaValues{ir}.t;
    medianAmps = tdpaValues{ir}.pMedianLog;
    ts         = tdpaValues{ir}.ts_arrival;
    
    subplot(ax.a1); hold on; grid on;
    plot(t,medianAmps,'lineStyle','-','color',colours{ir},'lineWidth',2);
end
xlabel('Time since P-onset [sec]','fontSize',ftSize)
set(gca,'fontSize',ftSize,'xlim',[tRange(1) tRange(2)])
set(gca,'ylim',[-9 -1],'xscale','log')
ylabel(['log_{10} (Pd) [m]'],'fontSize',ftSize)
%title(sprintf('%s %s on records with %i<=R<%ikm, SNR>=%i, %s records',ornt,band,rRange,fig.snrMin,titleAppendix))


ylm = get(gca,'ylim');
dy  = 0.1;
for ir = 1:nr
    ts = tdpaValues{ir}.ts_arrival;
    plot(ts,ylm(1)+(nr-ir+1)*dy,'.','color',colours{ir},'markerSize',10)
end