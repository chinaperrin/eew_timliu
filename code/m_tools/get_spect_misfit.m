function [m1,r1] = get_spect_misfit(zTarget,zTraining,snippet,nsim,bands,wt)

% Find magnitudes of most similar traces only wrt/ bands specified in
% <bands>, e.g. [1,2,8,9], or [1:9]

global o minimax

ftSize = 12;

rmin = minimax.rmin;
rmax = minimax.rmax;
mmin = minimax.mmin;
mmax = minimax.mmax;
amin = minimax.amin;
amax = minimax.amax;

m_ctlg = zTarget.m;
r_ctlg = zTarget.hypDist;

amax_obs      = zTarget.amax{1}(bands,snippet);
amax_training = cell2mat(cellfun(@(x) x(bands,snippet)', zTraining.amax,'UniformOutput',0));

nt         = size(amax_training,1);
misfit     = amax_training - repmat(amax_obs',nt,1);

% Cumulative misfit over all frequency bands
misfit_l2  = sum(misfit.^2,2);
%misfit_l1  = sum(abs(misfit),2);

% L2 norm cumulative misfit
[~,idx_std] = sort(misfit_l2);
mVal        = zTraining.m            (idx_std(1:nsim));
rVal        = log10(zTraining.hypDist(idx_std(1:nsim)));

% Unweighted parameter estimates
m0 = mean(mVal);
r0 = mean(rVal);

% Take best estimate to be weighted arithmetic mean
% Make sparse data points count more: replicate them according to weight
for ival = 1:nsim
    m       = mVal(ival);
    r       = rVal(ival);
    weight  = wt(ival);
    mVal    = [mVal; repmat(m,weight-1,1)];
    rVal    = [rVal; repmat(r,weight-1,1)];
end

m1    = mean(mVal);
r1log = mean(rVal);
r1    = 10^r1log;
% % L1 norm cumulative misfit 
% [~,idx_std_l1] = sort(misfit_l1);
% mVal_l1        = zTraining.m            (idx_std_l1(1:nsim));
% rVal_l1        = log10(zTraining.hypDist(idx_std_l1(1:nsim)));
% m2             = mean(mVal_l1);
% r2             = mean(rVal_l1);



if (o.plotSpectMisfit)
    
    figure(123); clf;
    whitebg(gcf,'k')
    set(gcf,'InvertHardcopy','off')

    % Read and plot waveform
    subplot(4,1,1); hold on; grid on
    [zS,meta] = read_any_trace(zTarget.fullName{1},zTarget,1);
    zPxIdx      = meta.ppxIdx;
    tz          = meta.t;
    ztpx        = tz(zPxIdx);
    plot(tz,zS.vel,'lineWidth',1)
    ymax = 1.1*max(abs(zS.vel));
    line([ztpx ztpx],[-ymax ymax],'Color','r','lineWidth',2)
    xlabel('Time [sec]','fontSize',ftSize)
    ylabel('Velocity [m/s]','fontSize',ftSize)
    set(gca,'fontSize',ftSize)
    set(gca,'ylim',[-ymax ymax])
    set(gca,'xlim',[min(tz) max(tz)])
    
    
    % Parameter estimation plot
    subplot(4,1,2:3); hold on
    %scatter(mVal,rVal,2,[1 1 1],'marker','.')
    scatter(10.^rVal,mVal)
    l1=line([rmin rmax],    [m_ctlg m_ctlg],'color','w');
       line([r_ctlg r_ctlg],[mmin mmax],    'color','w');
    p2=plot(10^r0,m0,'pw','markerSize',18,'markerFaceColor','w');
    p1=plot(10^r1,m1,'pw','markerSize',15,'markerFaceColor','r');
    lgnd=legend([p1;p2;l1],'Weighted L2-norm estimate','Unweighted L2-norm estimate','Catalog values');
    if (m_ctlg>5.5); set(lgnd,'location','southeast','fontSize',ftSize); end
    whitebg(gcf,'k'); set(gcf,'InvertHardcopy','off')
    %title('Magnitude and distance values of traces with most similar "spectra"','fontSize',ftSize)
    xlabel(['Epicentral Distance [km]'],'fontSize',ftSize)
    ylabel(['Magnitude'],'fontSize',ftSize)
    set(gca,'ylim',[mmin mmax],'xlim',[rmin rmax],'fontSize',ftSize)
    
    % Most similar SMA curves
    subplot(4,1,4); hold on; grid on
    p3=plot(amax_training(idx_std(1:nsim),:)');
       plot(amax_obs,'-k','lineWidth',5);
    p4=plot(amax_obs,'-r','lineWidth',2);
    %set(gca,'yscale','log','ylim',[10^amin 10^amax])
    %set(gca,'ytick',logspace(amin,amax,5))
    ylabel('Amplitudes [m/s]','fontSize',ftSize)
    legend([p3(1);p4],'most similar traces','target trace')
end
