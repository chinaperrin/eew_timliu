function [repickList] = run_manual_repick_loop(repickList,pxListName,px,opts)

% Manually go through repickList. If something, e.g. a pick needs to be
% changed, modify correspoding list entry, then run add2pxList.m. If you
% only add changes to repickList the changes will be overwritten when the
% TraceList is loaded the next time.


fprintf(1,'make sure you understand workflow for manual picking!\n')
% Case 0: Dont recompute Pd & Co. simply rerun picker and then add to
%                        pxlist. will be used in next wformProc iteration.
% Case 1: Wrong pick
%         option i)     - set tppx0
%                       - rerun picker --> ppxIdxNew, tppxNew
%                       - recompute acc/vel/dsp with ppxIdxNew
%                       - recompute peak amps with ppxIdxNew and newly
%                         processed acc/vel/dsp --> pd
%                       - compute ds
%                       - shift pd from two steps earlier
%
%       option ii)      - enter tppxNew manually
%                       - compute ppxIdxNew
%                       - proceed as before
%
% Case 2: only need to shift Pd
%                       - reprocess data to get acc/vel/dsp, using ppxIdx from list
%                       - compute ds
%                       - shift pd from repickList.pd{itr}


global fOrder fMode snpLength iN

fprintf(1,'Loading the latest tdpa_outMat-file. You may want to recompute it by running cf_peak_amps_1subplot.m\n')
%load('tmp/tdpa_outMat.mat')

fig.print = false;

ntap         = px.Param.ntap;
fUp_px       = px.Param.fUp_px;
fLow_px      = px.Param.fLow_px;
fLow_prefilt = px.Param.fLow_prefilt;

ntr = numel(repickList.m);

for itr = 1:ntr
    
    fprintf(1,sprintf('\n%i/%i\n',itr,ntr))
    fullName = repickList.fullName{itr};
    repickList.printSingleTraceSummary(itr)
    
    % Read trace
    [s,meta] = read_any_trace(fullName,repickList,opts);
    t        = meta.t;
    ppxIdx   = meta.ppxIdx;
    tppx     = t(ppxIdx);
    sr       = meta.sr;
    ppxIdx0  = ppxIdx - sr;
    isSM     = (strcmp(repickList.instrCode{itr},'L') ||strcmp(repickList.instrCode{itr},'N') ||strcmp(repickList.instrCode{itr},'G'));
    
    % Remove early mean, taper and pre-filter for picking ...
    sm   = s.raw - mean(s.raw(1:2*ntap));
    stap = taper(sm,sr,ntap);
    sb   = bworth(stap,sr,[fLow_prefilt,fUp_px],'band',2,'causal');       
    sh   = bworth(stap,sr,fLow_prefilt         ,'high',2,'causal');       

    % Look at waveform and existing pick
    repickList.printSingleTraceSummary(itr)
    clear sig allPx
    sig.sraw       = s.raw; % Raw, only gain corrected
    sig.sb         = sb;  % ... mean removed, tapered, highpass filtered
    sig.sh         = sh;   % ... mean removed, tapered, highpass filtered
    allPx.tppxList = tppx;
    plot_wform_components(sig,t,allPx,[tppx-.2 tppx+.05],[])
    
    % 3. Plot single record tdpa curvea
    %[h,ax] = prepare_singleRecTdpa_figure(outMat,[0 10]);
    %subplot(ax.a1); hold on; grid on;
    %plot_singleRec_tdpa(s,t,itr,repickList,tppx,ax,fig)
    
     
    
    % If something is not right ... xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    o.neverTrue = false; 
    if o.neverTrue
        
        % summarised in m_stats/tmp.m for easier manual running
        
        % ... reprocess signal for picking
        if isSM; [acc,vel,dsp] = sm2accVelDsp(sh,sr,ppxIdx,'afterPx',fLow_prefilt,2,fMode);
        else     [acc,vel,dsp] = bb2accVelDsp(sh,sr,ppxIdx,'afterPx',fLow_prefilt,2,fMode);
        end
        abpx = bworth(acc,sr,[fLow_px,fUp_px],'band',2,'causal');       % Bandpass filter for picking
        %alpx               = bworth(acc,sr,[fUp_px]        ,'low',2,'causal');       % Bandpass filter for picking

        % ... repick signal ...............................................
        %px.Opt.animatePx=1; px.Opt.plotPx=1;
        %px.Opt.animatePx=0; px.Opt.plotPx=0;
        tppx0                  = 10;
        [~,ppxIdx0]            = min(abs(t-tppx0));
        [ppxIdxNew,snrNew,~,~] = SBPx(abpx,1/sr,ppxIdx0,px.Param,px.Weight,px.Opt);
        tppxNew                = t(ppxIdxNew);
        if isempty(tppxNew); tppxNew = 9999; end
        fprintf(1,[' snr: ',num2str(snrNew,'%4.0f')])

        %tppxNew       = 4.78;
        [~,ppxIdxNew] = min(abs(t-tppxNew));

        % ... reprocess signal when there is a new pick
        if isSM; [acc,vel,dsp] = sm2accVelDsp(sh,sr,ppxIdxNew,'afterPx',fLow_prefilt,2,fMode);
        else     [acc,vel,dsp] = bb2accVelDsp(sh,sr,ppxIdxNew,'afterPx',fLow_prefilt,2,fMode);
        end
        abpx = bworth(acc,sr,[fLow_px,fUp_px],'band',2,'causal');       % Bandpass filter for picking

        
        clear sig allPx
        sig.sraw       = s.raw; % Raw, only gain corrected
        sig.sbpx       = abpx;  % ... mean removed, tapered, highpass filtered
        sig.dsp        = dsp;   % ... hp, then integrated, then bandpass filtered
        allPx.tppxList = tppx;
        allPx.tppxNew  = tppxNew;
        %plot_wform_components(sig,t,allPx,[tppx-1 tppx+.1],[])
        plot_wform_components(sig,t,allPx,[tppx-.2 tppx+.05],[])
        
        
        
        % ... recompute peak amps .........................................
        if repickList.m(itr)>7; tmax = 60; else tmax = 20; end
        [pa,pna,paIdx,~,~] = measure_tdpa(acc,ppxIdxNew,sr,opts.tnoise,opts.tmax,snpLength,[]);
        [pv,pnv,pvIdx,~,~] = measure_tdpa(vel,ppxIdxNew,sr,opts.tnoise,opts.tmax,snpLength,[]);
        [pd,pnd,pdIdx,~,~] = measure_tdpa(dsp,ppxIdxNew,sr,opts.tnoise,opts.tmax,snpLength,[]);
        pa                 = trim_amax(pa ,opts.nsnpmin);
        pv                 = trim_amax(pv ,opts.nsnpmin);
        pd                 = trim_amax(pd ,opts.nsnpmin);
        commentPx          ='manual pick';
        commentPx          ='automatically repicked';
        [snrNew,~,~,~]     = get_snr(abpx,ppxIdxNew,sr,px.Param.signalWindow,px.Param.gapWindow,px.Param.noiseWindow);

        %         % ... recompute pick delay ........................................
        %         copied to appendix
        
        % Save new values? ................................................
        repickList.pa{itr}       = single(pa);
        repickList.pv{itr}       = single(pv);
        repickList.pd{itr}       = single(pd);
        repickList.tppx(itr)     = tppxNew;
        repickList.ppxIdx(itr)   = ppxIdxNew;
        repickList.snr(itr)      = snrNew;
        repickList.comment{itr}  = commentPx;
        
        %         %repickList.pd{itr}       = single(pd_shifted);
        %         repickList.var2{itr}(1)  = ds;
        %         fprintf(1,'modify this statement after processing run is through\n')
                
        % SAVE NEW VALUES BY ADDING THEM TO PICKLIST
        pxList = add2pxList(itr,repickList,pxListName);
        
        
        % Alternatively: add trace to blackList
        blackListName = sprintf('blackList_i%i.mat',iN);
        comment     = 'no pre-pick signal';
        comment     = 'emergent signal';
        comment     = 'unclear pick';
        comment     = 'high noise';
        comment     = 'weird signal';
        comment     = 'in foreshock coda';
        [blackList] = add2blackList(itr,repickList,comment,blackListName);
    end
    1+1;
end



%% APPENDIX
%         % ... recompute pick delay ........................................
%         nrnd           = 100;
%         n              = 3;
%         n              = 2;
%         opts.plotSynth = 0;
%
%         abpx               = bworth(acc,sr,[fLow_px,fUp_px],'band',2,'causal');       % Bandpass filter for picking
%         [deltasVect,nAmps] = get_median_pick_delay(abpx,meta,n,1,nrnd,px,opts);
%         ds                 = round(median(deltasVect (deltasVect ~=99999)));
%         dsOld              = repickList.var2{itr}(1);
%         pd_shifted         = [zeros(1,ds), repickList.pd{itr}];
%         pd_shifted         = [zeros(1,ds), pd];
%
%         hf = figure(12); clf; box on; grid on; whitebg('w')
%         replot_tdpa_medians(outMat,[0 100]);
%         tsnp   = snpLength:snpLength:numel(pd_shifted)*snpLength;
%         plot(tsnp,log10(pd_shifted),'-xw','lineWidth',2)
%         plot(tsnp,log10(pd_shifted),'-xr','lineWidth',1)
%
%         plot_wform_and_tdpa(dsp,t,tppx-ds/sr   ,[],[],pd_shifted,snpLength,9001)
%         plot_wform_and_tdpa(dsp,t,tppxNew-ds/sr,[],[],pd_shifted,snpLength,9001)
