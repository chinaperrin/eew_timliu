function [TraceList] = overwritePx(TraceList,pxListFullName,px,opts)

global snpLength fOrder fMode

load(pxListFullName)


npx     = numel(pxList.m);
nsnpmin = 1000;            

% Unpack input arguments
ntap         = px.Param.ntap;
fUp_px       = px.Param.fUp_px;
fLow_px      = px.Param.fLow_px;
fLow_prefilt = px.Param.fLow_prefilt;
sWindow      = px.Param.signalWindow;
gWindow      = px.Param.gapWindow;
nWindow      = px.Param.noiseWindow;


for ipx = 1:npx
    
    % Find traceIdx and indices of other comps of same record
    fullName         = pxList.fullName{ipx};
    [idxZ,idxE,idxN] = find_corec_idx(fullName,TraceList,0);
    m                = pxList.m(ipx);
    tnoise           = numel(pxList.dspNoise{ipx})/pxList.sRate(ipx);
    sr               = pxList.sRate(ipx);
    isSM             = (strcmp(pxList.instrCode{ipx},'L') ||strcmp(pxList.instrCode{ipx},'N') ||strcmp(pxList.instrCode{ipx},'G'));
    
    
    % Warning
    if ~strcmp(pxList.orntCode{ipx},'Z')
        fprintf(1,'8UNG: this is a horizontal record. pxList should only contain vertical records.\n')
        pause(.5)
    end
    
    ntrfound = sum(cellfun(@(x) ~isempty(x), {idxZ; idxE; idxN}));
    if ntrfound==3
        
        r = TraceList.hypDist(idxZ);
        fprintf(1,sprintf('\n%i/%i: Overwriting picks in TraceList for following traces  (m%3.1f@%5.1fkm):',ipx,npx,m,r))
        TraceList.fullName([idxZ; idxE; idxN])
        
        % Compare new and old picks
        tppxOrg = TraceList.tppx([idxZ,idxE,idxN]);
        tppxNew = pxList.tppx(ipx);
        
        
        % Z: Read & process VERTICAL trace ............................
        [s,meta]      = read_any_trace(fullName,TraceList,opts);
        t             = meta.t;
        [~,ppxIdxZ]   = min(abs(tppxNew-t)); %ppxIdx   = meta.ppxIdx;
        sr            = meta.sr;
        sm            = s.raw - mean(s.raw(1:ppxIdxZ));
        stap          = taper(sm,sr,ntap);
        sh            = bworth(stap,sr,fLow_prefilt    ,'high',2,fMode);
        sb            = bworth(stap,sr,[fLow_px,fUp_px],'band',2,'causal');
        
        if isSM; [acc,vel,dsp] = sm2accVelDsp(sb,sr,ppxIdxZ,'afterPx',fLow_prefilt,2,fMode);
        else     [acc,vel,dsp] = bb2accVelDsp(sb,sr,ppxIdxZ,'afterPx',fLow_prefilt,2,fMode);
        end

        % Recompute peak amps
        if m<=7; tmax=20; else tmax = 60; end
        [pa,~,~,~,~] = measure_tdpa(acc,ppxIdxZ,sr,tnoise,tmax,snpLength,[]);
        [pv,~,~,~,~] = measure_tdpa(vel,ppxIdxZ,sr,tnoise,tmax,snpLength,[]);
        [pd,~,~,~,~] = measure_tdpa(dsp,ppxIdxZ,sr,tnoise,tmax,snpLength,[]);
        pa           = trim_amax(pa,nsnpmin);    % Trim time series when they no longer change
        pv           = trim_amax(pv,nsnpmin);
        pd           = trim_amax(pd,nsnpmin);
        [snr,~,~,~]  = get_snr(acc,ppxIdxZ,sr,sWindow,gWindow,nWindow);
        
        % Overwrite fields
        TraceList.ppxIdx     (idxZ) = ppxIdxZ;
        TraceList.tppx       (idxZ) = tppxNew;
        TraceList.snr        (idxZ) = snr;
        TraceList.pa         {idxZ} = single(pa);
        TraceList.pv         {idxZ} = single(pv);
        TraceList.pd         {idxZ} = single(pd);
        TraceList.dataSetName{idxZ} = sprintf('%s_manPx',TraceList.dataSetName{idxZ});
        
        % verify that new px is better (visually)
        if opts.plotNewPx
            
            pxList.printSingleTraceSummary(ipx)
            
            sig.sraw      = s.raw;  % Raw, only gain corrected
            sig.sb        = sb;     % ... dito but bandpass filtered
            sig.sh        = sh;     % ... dito but highpass filtered
            allPx.tppxOrg = tppxOrg(1);
            allPx.tppxNew = tppxNew;
            plot_wform_components(sig,t,allPx,[tppxNew-.3 tppxNew+.5],[])
            1+1;
            
            % Manually repick if necessary
            o_neverTrue = false;
            if o_neverTrue
                
                % Find new pick
                tppx0                  = 20.6-1;
                [~,ppxIdx0]            = min(abs(t-tppx0));
                [ppxIdxNew,snrNew,~,~] = SBPx(sb,1/sr,ppxIdx0,px.Param,px.Weight,px.Opt);
                tppxNew                = t(ppxIdxNew);
                
                % Recompute peak amps
                if m<=7; tmax=20; else tmax = 60; end
                [pa,~,~,~,~] = measure_tdpa(acc,ppxIdxNew,sr,tnoise,tmax,snpLength,[]);
                [pv,~,~,~,~] = measure_tdpa(vel,ppxIdxNew,sr,tnoise,tmax,snpLength,[]);
                [pd,~,~,~,~] = measure_tdpa(dsp,ppxIdxNew,sr,tnoise,tmax,snpLength,[]);
                pa           = trim_amax(pa,nsnpmin);    % Trim time series when they no longer change
                pv           = trim_amax(pv,nsnpmin);
                pd           = trim_amax(pd,nsnpmin);
                [snr,~,~,~]  = get_snr(acc,ppxIdxNew,sr,sWindow,gWindow,nWindow);
                
                % Overwrite fields in TraceList
                TraceList.ppxIdx     (idxZ) = ppxIdxNew;
                TraceList.tppx       (idxZ) = tppxNew;
                TraceList.snr        (idxZ) = snr;
                TraceList.pa         {idxZ} = single(pa);
                TraceList.pv         {idxZ} = single(pv);
                TraceList.pd         {idxZ} = single(pd);
                TraceList.dataSetName{idxZ} = sprintf('%s_manPx',TraceList.dataSetName{idxZ});

                % Overwrite fields in pxList
                pxList.ppxIdx     (ipx) = ppxIdxNew;
                pxList.tppx       (ipx) = tppxNew;
                pxList.snr        (ipx) = snr;
                pxList.pa         {ipx} = single(pa);
                pxList.pv         {ipx} = single(pv);
                pxList.pd         {ipx} = single(pd);
                pxList.dataSetName{ipx} = sprintf('%s_manPx',TraceList.dataSetName{idxZ});
            
                save(pxListFullName,'pxList')
            end 
        end
        
        
        % E: Read & process EAST trace ................................
        [s,meta]      = read_any_trace(TraceList.fullName{idxE},TraceList,opts);
        t             = meta.t;
        [~,ppxIdxE]   = min(abs(tppxNew-t));    % Pick index on the time vector of the actual
        sm            = s.raw - mean(s.raw(1:ppxIdxE));
        stap          = taper(sm,sr,ntap);
        sb            = bworth(stap,sr,[fLow_px,30],'band',2,'causal');
        [acc,vel,dsp] = sm2accVelDsp(sb,sr,ppxIdxE,'afterPx',fLow_prefilt,2,fMode);
        
        % Recompute peak amps
        if m<=7; tmax=20; else tmax = 60; end
        [pa,~,~,~,~] = measure_tdpa(acc,ppxIdxE,sr,tnoise,tmax,snpLength,[]);
        [pv,~,~,~,~] = measure_tdpa(vel,ppxIdxE,sr,tnoise,tmax,snpLength,[]);
        [pd,~,~,~,~] = measure_tdpa(dsp,ppxIdxE,sr,tnoise,tmax,snpLength,[]);
        pa           = trim_amax(pa,nsnpmin);    % Trim time series when they no longer change
        pv           = trim_amax(pv,nsnpmin);
        pd           = trim_amax(pd,nsnpmin);
        [snr,~,~,~]  = get_snr(acc,ppxIdxE,sr,sWindow,gWindow,nWindow);
         
        % Overwrite fields
        TraceList.ppxIdx     (idxE) = ppxIdxE;
        TraceList.tppx       (idxE) = tppxNew;
        TraceList.snr        (idxE) = snr;
        TraceList.pa         {idxE} = single(pa);
        TraceList.pv         {idxE} = single(pv);
        TraceList.pd         {idxE} = single(pd);
        TraceList.dataSetName{idxE} = sprintf('%s_manPx',TraceList.dataSetName{idxE});
        
         
        
        % N: Read & process NORTH trace ...............................
        [s,meta]      = read_any_trace(TraceList.fullName{idxN},TraceList,opts);
        t             = meta.t;
        [~,ppxIdxN]   = min(abs(tppxNew-t));    % Pick index on the time vector of the actual
        sm            = s.raw - mean(s.raw(1:ppxIdxN));
        stap          = taper(sm,sr,ntap);
        sb            = bworth(stap,sr,[fLow_px,30],'band',2,'causal');
        [acc,vel,dsp] = sm2accVelDsp(sb,sr,ppxIdxN,'afterPx',fLow_prefilt,2,fMode);
        
        % Recompute peak amps
        if m<=7; tmax=20; else tmax = 60; end
        [pa,~,~,~,~] = measure_tdpa(acc,ppxIdxN,sr,tnoise,tmax,snpLength,[]);
        [pv,~,~,~,~] = measure_tdpa(vel,ppxIdxN,sr,tnoise,tmax,snpLength,[]);
        [pd,~,~,~,~] = measure_tdpa(dsp,ppxIdxN,sr,tnoise,tmax,snpLength,[]);
        pa           = trim_amax(pa,nsnpmin);    % Trim time series when they no longer change
        pv           = trim_amax(pv,nsnpmin);
        pd           = trim_amax(pd,nsnpmin);
        [snr,~,~,~]  = get_snr(acc,ppxIdxN,sr,sWindow,gWindow,nWindow);
        
        % Overwrite fields
        TraceList.ppxIdx     (idxN) = ppxIdxN;
        TraceList.tppx       (idxN) = tppxNew;
        TraceList.snr        (idxN) = snr;
        TraceList.pa         {idxN} = single(pa);
        TraceList.pv         {idxN} = single(pv);
        TraceList.pd         {idxN} = single(pd);
        TraceList.dataSetName{idxN} = sprintf('%s_manPx',TraceList.dataSetName{idxN});
        
    else
        fprintf(1,'At least one component of pxList entry not found in traceList. \nNot overwriting corresponding fields in TraceList.\n')
    end
end
