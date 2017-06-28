function trList = recompute_Pd_with_added_noise(trList,params,opt)

global iN fOrder fMode snpLength dataSetNames

%dbstop if error



% Should be passed on as global variables
%snpLength = trList.prop.snpLength
%fMode     = 'causal'


% PICKER (see SBPx.m for meaning of parameters)
px           = load_picker_settings('production_short');
fLow_prefilt = px.Param.fLow_prefilt;
fLow_px      = px.Param.fLow_px;
fUp_px       = px.Param.fUp_px;
ntap         = px.Param.ntap;

nrnd         = params.nrnd;
nzcMin       = params.nzcMin;
tmax         = params.tmax; 
tmax_bak     = tmax;
nsnpmin      = params.nsnpmin;
nsnpmin_bak  = nsnpmin;

outListFullName = sprintf('/scratch/memeier/fbout/i%i/zLists/noiseScaled/new/scaledNoiseList.mat',iN);
%outVarFullName  = sprintf('/scratch/memeier/fbout/i%i/zLists/noiseScaled/new/scaledNoiseList_vars.mat',iN);

% Read trace
ntr         = numel(trList.m);
wforms      = cell(ntr,1);
scaleFactor = zeros(ntr,1);

for itr = 1:ntr

    %if (ismember(itr,hundreds) |itr==1 |itr==ntr); print_iteration_numbers(itr,ntr,'hundreds'); end
    print_iteration_numbers(itr,ntr,'tens')
    if opt.debug; trList.printSingleTraceSummary(itr); end
        
    
    [S,meta]  = read_any_trace(trList.fullName{itr},trList,opt);
    m         = trList.m(itr);
    t         = meta.t;
    sr        = meta.sr;
    ppxIdx    = trList.ppxIdx(itr);          % Read original, uncorrected pick
    tppx      = t(ppxIdx);
    instrCode = trList.instrCode(itr);
    isSM      = (strcmp(instrCode,'L') ||strcmp(instrCode,'N') ||strcmp(instrCode,'G'));
    ns        = numel(S.raw);
    
    % Prepare vertical trace for picking
    sm   = S.raw - mean(S.raw(1:2*ntap));
    stap = taper(sm,sr,ntap);
    stb  = bworth(stap,sr,[fLow_prefilt,30],'band',2,fMode);
    
    if ~isSM; [acc,~,~] = bb2accVelDsp(stb,sr,ppxIdx,'allWform',fLow_prefilt,2,fMode);
    else      acc       = stb;
    end
    
    if opt.debug; h = plot_wform(S.raw,t,tppx,[],[tppx-2 tppx+24],20); end
    
    % Cut out 'reliable' noise part
    if     ppxIdx>2000; isn = 1000;          % If there is a lot of pre-signal waveform, cut out more...
    elseif ppxIdx>1000; isn = 500;
    else                isn = 2*ntap;
    end
    ien   = round(.95*ppxIdx);
    snRaw = acc(isn:ien);
    snRaw = snRaw - mean(snRaw);
    %snRaw = bworth(snRaw,sr,fLow_px,'high',2,'causal');
    if opt.debug; clf; hold on; plot(acc,'y'); plot(isn:ien,snRaw,'r'); set(gca,'xlim',[1,ppxIdx+40]); end
    
    % Concatenate replications of snRaw: i) find last three zero crossings,
    % ii) find last positive maximum in between these, iii) do same at
    % beginning of snRaw, iv) use the segment between those indices 
    % (=snSegment) for concatenation.
    if opt.debug; oc.plot = true; else oc.plot=false; end
    crIdx = get_zero_crossings(snRaw,oc);
    nzc   = numel(crIdx);
    
    if nzc>=nzcMin
        [~,iFirstMax]   = max(snRaw(1:crIdx(2)+1));
        [~,iLastMaxRel] = max(snRaw(crIdx(end-2):crIdx(end)));
        iLastMax        = crIdx(end-2)+iLastMaxRel-1;
        snSegment       = snRaw(iFirstMax:iLastMax);
        nsseg           = numel(snSegment);
        
        % Replicate snSegment to be long enough for random onset index sampling
        nspx        = ppxIdx+5*sr;   % Length of signal to which noise is added
        if nspx>ns; nspx = ns; end
        nscatTarget = 5*nspx;        % Make length of noise signal much longer
        nrep        = ceil(nscatTarget/nsseg);
        snCat       = [snRaw(1:iLastMax); repmat(snSegment,nrep,1)];
        idxCat      = [iLastMax;          iLastMax+(1:nrep)'*nsseg];
        nscat       = numel(snCat);
        if opt.debug; clf; hold on; plot(snCat,'y'); plot(snRaw,'r'); plot(idxCat,0,'dm','markerSize',12); end
        
        % Scale concatenated noise vector snCat to target amplitude
        p95      = prctile(snCat,95);
        p95list  = trList.accNoise{itr}(3);
        scFactor = params.p95target/p95;
        if scFactor<1; 
            scFactor = 0; 
        end    % if noise level is higher than target, don't add any noise.
        snScaled = scFactor*snCat;
        if opt.debug; gcf; hold on; plot(snScaled,'w'); end
        
        % Choose random starting point and add to signal
        ppxIdxMat = zeros(nrnd,1);
        
        for irnd = 1:nrnd
            
            imax   = nscat-nspx-1;
            if imax<1; imax=1; end
            isrand = randi(imax,1,1);
            ierand = isrand+nspx-1;
            if ierand>nscat; ierand=nscat; end
            noise     = snScaled(isrand:ierand);
            accOrg    = acc(1:nspx);
            accNoised = accOrg+noise;
            p95check  = prctile(accNoised(ntap:round(.95*ppxIdx)),95);         % Check p95 of noise after summation
            
            if opt.debug;
                clf; hold on; plot(accOrg,'r'); plot(noise,'y'); plot(accNoised,'m')
                fprintf(1,sprintf('p95: %em/s/s before, %em/s/s after scaling, target=%em/s/s\n',p95,p95check,params.p95target))
            end
            
            %px.Opt.plotPx=1
            apx              = bworth(accNoised,sr,[fLow_px fUp_px],'band',2,'causal');
            [ppxIdx_n,~,~,~] = SBPx(apx,1/sr,ppxIdx-sr,px.Param,px.Weight,px.Opt);
            dpxIdxMax = .5*sr;
            dpxIdx    = abs(ppxIdx_n-ppxIdx);
            flg_pxFound  = ~isempty(ppxIdx_n);
            if flg_pxFound; flg_pxTooOff = (dpxIdx>dpxIdxMax);
                if ~flg_pxTooOff; ppxIdxMat(irnd) = ppxIdx_n;
                else              ppxIdxMat(irnd) = 99999;
                end
            else            ppxIdxMat(irnd)  = 99999;
            end
        end
        
        % Save last waveform and noise scaling factor
        wforms{itr}      = accNoised;
        scaleFactor(itr) = scFactor;
        
        medPpxIdx = round(median(ppxIdxMat(ppxIdxMat~=99999)));
        minPxIdx  = min(ppxIdxMat(ppxIdxMat~=99999));
        trList.var1{itr}.medPpxIdx = medPpxIdx;
        trList.var1{itr}.ppxIdxMat = ppxIdxMat;
    else
        fprintf(1,'Not enough zero crossings...\n')
    end
    
    if ( ~isnan(medPpxIdx) &&(nzc>=nzcMin) &&numel(snRaw)>4*ntap )
        
        % Recompute Pd-curves: Remove early mean, taper and pre-filter waveform (2-pole high pass) ...
        sm   = S.raw - mean(S.raw(1:medPpxIdx));
        stap = taper(sm,sr,ntap);
        sb   = bworth(stap,sr,[fLow_prefilt,30],'band',2,fMode);
        
        intPxIdx = medPpxIdx-params.nsprepx;
        if intPxIdx<1; intPxIdx=1; end
        if isSM; [acc2,~,dsp2] = sm2accVelDsp(sb,sr,intPxIdx               ,'afterPx',fLow_prefilt,2,fMode);
                 [~   ,~,dsp3] = sm2accVelDsp(sb,sr,minPxIdx-params.nsprepx,'afterPx',fLow_prefilt,2,fMode);
        else     [acc2,~,dsp2] = bb2accVelDsp(sb,sr,intPxIdx               ,'afterPx',fLow_prefilt,2,fMode);
                 [~   ,~,dsp3] = bb2accVelDsp(sb,sr,minPxIdx-params.nsprepx,'afterPx',fLow_prefilt,2,fMode);
        end
        %[acc3,~,dsp3] = sm2accVelDsp(sb,sr,medPpxIdx-100,'afterPx',fLow_prefilt,2,fMode);
        
        % Recompute pd-curve
        if m>7; params.tmax=60; nsnpmin=1e4; end
        [pd ,~,~,~,~]    = measure_tdpa(dsp2,medPpxIdx,sr,params.tnoise,params.tmax,snpLength,[]);
        [pd0,~,~,~,~]    = measure_tdpa(dsp3,minPxIdx ,sr,params.tnoise,params.tmax,snpLength,[]);
        [pa ,~,~,~,~]    = measure_tdpa(acc2,medPpxIdx,sr,params.tnoise,params.tmax,snpLength,[]);
        pd               = trim_amax(pd ,nsnpmin);
        pd0              = trim_amax(pd0,nsnpmin);
        pa               = trim_amax(pa ,nsnpmin);
        trList.var2{itr} = single(pd );
        trList.var3{itr} = single(pa );
        trList.var4{itr} = single(pd0);
            
        % Compute pd-curve on noise signal
        noisePpxIdx          = ntap;
        stapN                = taper(snRaw,sr,ntap);
        sbN                  = bworth(stapN,sr,[fLow_prefilt,30],'band',2,fMode);
        [accN,~,dspN]        = sm2accVelDsp(sbN,sr,noisePpxIdx,'afterPx',fLow_prefilt,2,fMode);
        [pdN,~,~,~,~]        = measure_tdpa(dspN,noisePpxIdx,sr,params.tnoise,params.tmax,snpLength,[]);
        [paN,~,~,~,~]        = measure_tdpa(accN,noisePpxIdx,sr,params.tnoise,params.tmax,snpLength,[]);
        pdN                  = trim_amax(pdN,nsnpmin);
        paN                  = trim_amax(paN,nsnpmin);
        trList.var5{itr}.pdN = single(pdN);
        trList.var5{itr}.paN = single(paN);
        
        % Save displacement after pick, interpolated to snpLength intervals
        tsnp                    = (1:nsnpmin)'*snpLength;
        dspSnp                  = interp1(t(medPpxIdx:end),dsp2(medPpxIdx:end),tsnp);
        trList.var5{itr}.dspSnp = single(dspSnp);
        % clf; plot(t,dsp2,'xr'); hold on; plot(tsnp,dspSnp,'-dk')
        
        %         isDsp = medPpxIdx;
        %         ieDsp = medPpxIdx+nsnpmin;
        %         if ieDsp>numel(dsp2); ieDsp=numel(dsp2); end
        %         trList.var5{itr}.absDsp = single(dsp2(isDsp:ieDsp));
        if m>7; params.tmax=tmax_bak; nsnpmin=nsnpmin_bak; end
    else
        trList.var2{itr} = nan;
    end
end


if opt.saveOut
    fprintf(1,'\nSaving out-file ... ')
    trList.prop.recompPd.params       = params;
    trList.prop.recompPd.opt          = opt;
    trList.prop.recompPd.px           = px;
    trList.prop.recompPd.wforms       = wforms;
    trList.prop.recompPd.scaleFactors = scaleFactor;
    trList.prop.recompPd.dataSetNames = dataSetNames;
    %save(outVarFullName ,'out')
    save(outListFullName,'trList')
    fprintf(1,'done\n ')
end