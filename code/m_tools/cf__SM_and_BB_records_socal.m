function cf__SM_and_BB_records_socal(catList)

% Use colocated BB and SM records to check the effect that two integration
% methods have on the Pd(t) curves.

opt_readWform.plotSynth  = 0;
opt_readWform.process    = 1;
opt_readWform.intMode    = 'afterPx';
opt_readWform.scp_wforms = true;

ftSize = 15;

fprintf(1,'Finding double records  ...')
matches = find_colocated_SM_BB_records(catList);
% matches_bak = matches;

% Remove multiple entries
lgc1          = (cellfun(@(x) numel(x), matches.bbIdx)>1);
lgc2          = (cellfun(@(x) numel(x), matches.smIdx)>1);
lgc           = logical(lgc1 |lgc2);
matches.bbIdx = matches.bbIdx(~lgc);
matches.smIdx = matches.smIdx(~lgc);

% Throw out all records/matches with hypDist>25km
% lgc           = catList.hypDist(cell2mat(matches.bbIdx))<=25;
% matches.bbIdx = matches.bbIdx(lgc);
% matches.smIdx = matches.smIdx(lgc);


%% Compare integration methods
nmatch = numel(matches.bbIdx);
for imatch = 1:nmatch
    
    fprintf(1,sprintf('%i/%i\n',imatch,nmatch))
    
    ibb = matches.bbIdx{imatch};
    ism = matches.smIdx{imatch};
    
    catList.printSingleTraceSummary(ibb)
    catList.printSingleTraceSummary(ism)
    %catList.printSingleTraceSummary(ism(2))
    
    % BROAD BAND
    [sbb,meta_bb] = read_any_trace(catList.fullName{ibb},catList,opt_readWform);
    t             = meta_bb.t;
    ppxIdx        = meta_bb.ppxIdx;
    tppx          = t(ppxIdx);
    
    % STRONG MOTION
    [ssm,meta_sm] = read_any_trace(catList.fullName{ism},catList,opt_readWform);
    t2            = meta_sm.t;
    ppxIdx2       = meta_sm.ppxIdx;
    tppx2         = t2(ppxIdx2);
    
    figure(11); clf; hold on; grid on; box on;
    pbb = plot(t ,sbb.acc,'b','lineWidth',2);
    psm = plot(t2,ssm.acc,'k','lineWidth',1);
    plot(tppx,0,'xr','lineWidth',2,'markerSize',12)
    set(gca,'xlim',[tppx-.2 tppx+.1])
    
    
    l1 = legend([pbb,psm],'BB record','SM record');
    set(l1,'fontSize',ftSize)
    1+1; 
end




% Plot curves
nmatch  = numel(matches.bbIdx);
for imatch = 1:nmatch
    
    ibb = matches.bbIdx{imatch};
    ism = matches.smIdx{imatch};
    
    catList.printSingleTraceSummary(ibb)
    catList.printSingleTraceSummary(ism)
    %catList.printSingleTraceSummary(ism(2))
    
    % Broad band
    [sbb,meta_bb]       = read_any_trace_bak(catList.fullName{ibb},catList,1);
    tbb                 = meta_bb.t;
    tppx_bb             = tbb(meta_bb.ppxIdx);
    sr_bb               = meta_bb.sr;
    ns                  = round(snpLength*meta_bb.sr);
    [pd_bb,~,pd_idx_bb] = measure_peakAmps_in_expanding_tWindows(s_bb.dsp,meta_bb.ppxIdx,sr_bb,0,100,ns);
    [pv_bb,~,pv_idx_bb] = measure_peakAmps_in_expanding_tWindows(s_bb.vel,meta_bb.ppxIdx,sr_bb,0,100,ns);
    
    % Strong motion
    [s_sm,meta_sm]      = read_any_trace_(catList.fullName{ism},catList,1);
    tsm                 = meta_sm.t;
    tppx_sm             = tsm(meta_sm.ppxIdx);
    sr_sm               = meta_sm.sr;
    ns                  = round(snpLength*meta_sm.sr);
    [pd_sm,~,pd_idx_sm] = measure_peakAmps_in_expanding_tWindows(s_sm.dsp,meta_sm.ppxIdx,sr_sm,0,100,ns);
    [pv_sm,~,pv_idx_sm] = measure_peakAmps_in_expanding_tWindows(s_sm.vel,meta_sm.ppxIdx,sr_sm,0,100,ns);
    
    % Asses the effect of highpass filtering i) after 1.int. & after 1. and
    % 2. int.
    ahih                  = bworth(s_sm.vel,meta_sm.sr,0.075,'high',fOrder,fMode);
    ahihi                 = acc2vel(ahih,meta_sm.sr,2);
    ahihih                = bworth(ahihi,meta_sm.sr,0.075,'high',fOrder,fMode);
    [pd_sm1,~,pd_idx_sm1] = measure_peakAmps_in_expanding_tWindows(ahihi,    meta_sm.ppxIdx,sr_sm,0,100,ns);
    [pd_sm2,~,pd_idx_sm2] = measure_peakAmps_in_expanding_tWindows(ahihih,meta_sm.ppxIdx,sr_sm,0,100,ns);
    
    [vh]                  = bworth(s_bb.vel,meta_bb.sr,0.075,'high',fOrder,fMode);
    vhi                   = acc2vel(vh,meta_bb.sr,2);
    vhih                  = bworth(vhi,meta_bb.sr,0.075,'high',fOrder,fMode);
    
    
    % Plot DISPLACEMENT
    figure(5); clf; hold on; grid on; box on;
    set(gca,'xlim',[tppx_bb-1, tppx_bb+1])
    
    p1 = plot(tbb,s_bb.dsp,'-k','lineWidth',2);
    plot(tppx_bb,0,'xb','lineWidth',2,'markerSize',13)
    plot(tbb(pd_idx_bb), pd_bb,':k')
    plot(tbb(pd_idx_bb),-pd_bb,':k')
    
    p2 = plot(tsm,s_sm.dsp,'-r','lineWidth',2);
    % plot(tsm,ahihi,'-b','lineWidth',1);
    plot(tppx_sm,0,'db','lineWidth',2,'markerSize',13)
    plot(tsm(pd_idx_sm), pd_sm,':r')
    plot(tsm(pd_idx_sm),-pd_sm,':r')
    
    l1 = legend([p1;p2],'Displacement (integrated BB trace)','Displacement (integrated SM trace)');
    set(l1,'fontSize',ftSize,'location','northWest')
    % print('-depsc2',givemeaname)
    
    
    % Plot VELOCITY
    figure(6); clf; hold on; grid on; box on;
    set(gca,'xlim',[tppx_bb-1, tppx_bb+1])
    
    p1b = plot(tbb,s_bb.vel,'-k','lineWidth',2);
    plot(tppx_bb,0,'xb','lineWidth',2,'markerSize',13)
    plot(tbb(pv_idx_bb), pv_bb,':k')
    plot(tbb(pv_idx_bb),-pv_bb,':k')
    
    p2b = plot(tsm,s_sm.vel,'-r','lineWidth',2);
    % plot(tsm,ahih,'-b','lineWidth',1);
    plot(tppx_sm,0,'db','lineWidth',2,'markerSize',13)
    plot(tsm(pv_idx_sm), pv_sm,':r')
    plot(tsm(pv_idx_sm),-pv_sm,':r')
    
    l2 = legend([p1b;p2b],'Velocity (BB trace)','Velocity (integrated SM trace)');
    set(l2,'fontSize',ftSize,'location','northWest')
    % print('-depsc2',givemeaname)
    
    
    % Plot pd curves on summary plot
    figure(hf);
    trel_bb = tbb(pd_idx_bb) - tppx_bb;
    trel_sm = tsm(pd_idx_sm) - tppx_sm;
    
    delete(px1)
    delete(px2)
    px1      = plot(trel_bb, log10(pd_bb),':xk','lineWidth',2);
    px2      = plot(trel_sm, log10(pd_sm),':xr','lineWidth',2);
    
    % Plot pd curves on separate plot
    figure(8); clf; hold on; grid on; box on;
    plot(trel_bb, log10(pd_bb),':k','lineWidth',4);
    plot(trel_sm, log10(pd_sm),':r','lineWidth',4);
    set(gca,'xscale','log','xlim',[8e-3 8],'ylim',[-7 -1])
    
    pause
    
    % Plot integrated SM record that has been filtered after first integration
    figure(5); hold on; plot(tsm,ahihi,'-m','lineWidth',1);                % --> sm trace
    figure(6); hold on; plot(tsm,ahih ,'-m','lineWidth',1);
    figure(8); hold on; plot(trel_sm, log10(pd_sm1),':m','lineWidth',4);
    
    % Plot integrated SM record that has been filtered after first integration
    pause
    
    % Plot integrated SM record that has been filtered after first and second integration
    figure(5); hold on; plot(tsm,ahihih,'-y','lineWidth',1);
    figure(5); hold on; plot(tbb,vhih  ,'-g','lineWidth',1);
    figure(8); hold on; plot(trel_sm, log10(pd_sm2),':y','lineWidth',4);
    pause
end