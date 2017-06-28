addpath(genpath('../'))

clear all



% Import traceList
opt.plotSynth  = 0;
opt.process    = 1;
opt.scp_wforms = 1;
opt.intMode    = 'afterPx';

iN             = 35;
fMode          = 'causal';
fOrder         = 4;
outDirAppendix = '_0p01';

dataSetBaseNames         = {'/scratch/memeier/data/japan/k_kik/M5p_z25km/out/'; ...
                            '/scratch/memeier/data/socal/scsn_120101_151231_M2_2p5/out/'};
outDirName               = sprintf('%s%s',fMode,outDirAppendix); 
outDirFullName           = sprintf('/scratch/memeier/fbout/i%i/%s/',35,outDirName);
dataSetNames             = strrep(dataSetBaseNames,'/out/',sprintf('/out/i%i/%s/',35,outDirName));
[TraceList,fc,snpLength] = import_trLists(dataSetNames,'traceList');



% Choose and read single waveform
itr = 120;
itr = 8100;

for itr = 1:100
    [S,meta]  = read_any_trace(TraceList.fullName{itr},TraceList,opt);
    m         = TraceList.m(itr);
    t         = meta.t;
    sr        = meta.sr;
    ppxIdx    = TraceList.ppxIdx(itr);          % Read original, uncorrected pick
    tppx      = t(ppxIdx);
    instrCode = TraceList.instrCode(itr);
    isSM      = (strcmp(instrCode,'L') ||strcmp(instrCode,'N') ||strcmp(instrCode,'G'));
    
    
    figure(22); clf; hold on;
    te   = tppx+.8;
    ts   = tppx-.1;
    
    ds          = TraceList.var2{itr}(1);
    ppxIdx_corr = ppxIdx-ds;
    tppx_corr   = t(ppxIdx_corr);
    plot(t,-S.dsp,'k');
    set(gca,'xlim',[ts te]);
    plot(tppx     ,0,'dk','markerFaceColor','y','markerSize',13)
    plot(tppx_corr,0,'dk','markerFaceColor','r','markerSize',13)
    
    pd   = TraceList.pd{itr};
    npd  = numel(pd);
    tpd  = tppx             +(1:npd)*snpLength;
    tpd2 = tppx-ds*snpLength+(1:npd)*snpLength;
    plot(tpd2,pd,'ok')
    
    % Reprocess waveform
    ntap = 10;
    fLow_prefilt = 0.075;
    sm   = S.raw - mean(S.raw(1:ppxIdx));
    stap = taper(sm,sr,ntap);
    ah   = bworth(stap,sr,fLow_prefilt,'high',2,fMode);
    
    ab               = bworth(stap,sr,[fLow_prefilt,30],'band',2,fMode);
    [acc2,vel2,dsp2] = sm2accVelDsp(ab,sr,ppxIdx_corr,'afterPx',fLow_prefilt,2,fMode);
    plot(t,abs(dsp2))
    set(gca,'ylim',[0 1e-5])
    
    % Find index of point where pd-cuve stays flat for at least <nflat> samples
    % for the first time
    nflat      = 10;
    dpd        = diff(pd);
    dpd        = [dpd(1),dpd];
    dpd_runsum = filter(ones(1,nflat),1,dpd);
    idx_flat   = find(dpd_runsum==0,1,'first')-nflat;
    idx_flat   = find(dpd_runsum==0);
    idx_flat   = idx_flat(idx_flat>nflat)
    idx_flat   = idx_flat(1)-nflat;
    plot(tpd(idx_flat),pd(idx_flat),'dk','markerFaceColor','r','markerSize',12)

    1+1;
end














figure(22); clf; hold on; 
subplot(3,1,1); hold on; plot(t,S.dsp,'k'); set(gca,'xlim',[ts te]);plot(tppx,0,'dk','markerFaceColor','y','markerSize',13) 
subplot(3,1,2); hold on; plot(t,S.vel,'k'); set(gca,'xlim',[ts te]);plot(tppx,0,'dk','markerFaceColor','y','markerSize',13)
subplot(3,1,3); hold on; plot(t,S.acc,'k'); set(gca,'xlim',[ts te]);plot(tppx,0,'dk','markerFaceColor','y','markerSize',13)
plot(t,abs(S.dsp),'r')
set(gca,'xlim',[ts te])
plot(tppx,0,'dk','markerFaceColor','y','markerSize',13)

ds = TraceList.var2{itr}(2);
dpd = diff(pd);
dpd = [dpd(1),dpd];



% Find index of point where pd-cuve stays flat for at least <nflat> samples
% for the first time
nflat=100;
dpd_runsum = filter(ones(1,nflat),1,dpd);
idx_flat   = find(dpd_runsum==0,1,'first')-nflat;
clf; hold on;
plot(tpd,pd)
plot(tpd(idx_flat),pd(idx_flat),'dk','markerFaceColor','r','markerSize',12)





figure(23); clf; hold on;
plot(t,S.acc,'k')
plot(t,abs(S.acc),'r')
set(gca,'xlim',[ts te])
plot(tppx,0,'dk','markerFaceColor','y','markerSize',13)