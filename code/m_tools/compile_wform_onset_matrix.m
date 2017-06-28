function [out] = compile_wform_onset_matrix(trList,params)

    fprintf(1,'8ung: Some hard-coded parameters at beginning of function. Check.\n')
    
    o.process   = 1;
    o.intMode   = 'afterPx';
    o.scpWforms = true;
    
    dtInp    = params.dtInp;             % Target sampling interval
    tSignal  = params.tSignal;    % Time in [s] after tppx that is saved to array
    tNoise   = params.tNoise;
    nsamples = params.nsamples;  % Max no of records in each m-bin
    
    mVect    = 2:1:8;
    
    sampleList = traceList(0);
    for im  = 1:numel(mVect)-1
        mLo = mVect(im);
        mUp = mLo+1;
        idx = find(trList.m>=mLo &trList.m<mUp);
        if numel(idx)>nsamples
            idx = datasample(idx,nsamples,'Replace',false);
        end
        tmpList = trList.selectSubList(idx);
        sampleList.appendList(tmpList);
    end
    
    
    %     idx80      = find(sampleList.sRate==80);
    %     sr80List   = sampleList.selectSubList(idx80);
    %     sampleList = sr80List;    
    %
    %     idx200     = find(sampleList.sRate==200);
    %     sr200List  = sampleList.selectSubList(idx200);
    %     sampleList = sr200List;
    
    ntr         = numel(sampleList.m);
    RAW         = cell (ntr,1);
    ACC         = cell (ntr,1);
    VEL         = cell (ntr,1);
    DSP         = cell (ntr,1);
    TIME        = cell (ntr,1);
    META.m      = zeros(ntr,1);
    META.r      = zeros(ntr,1);
    META.sr     = zeros(ntr,1);
    META.ppxIdx = cell (ntr,1);
    META.tppx   = zeros(ntr,1);
    
    
    for itr = 1:ntr
        
        
        % Read/process waveform, read meta-info
        [S,meta]  = read_any_trace(sampleList.fullName{itr},sampleList,o);
        tOrg      = meta.t;
        ppxIdxOrg = sampleList.ppxIdx(itr);
        tppxOrg   = tOrg(ppxIdxOrg);
        
        % Use modified pick onset! 
        
        fprintf(1,sprintf('%i/%i: sr = %i\n',itr,ntr,meta.sr))

        % Interpolate time series data
        dtOrg   = 1/meta.sr;
        srRatio = dtOrg/dtInp;        % Ratio between original and target sampling interval

        % Create interpolation time vector tInp
        diLo  = tNoise/dtInp;
        diUp  = tSignal/dtInp;
        tInp  = tppxOrg-diLo*dtInp:dtInp:tppxOrg+diUp*dtInp;
        nsInp = numel(tInp);
        
        % Pick index, pick time and pick delay shift on interpolated vector
        ppxIdxInp = diLo+1;
        tppxInp   = tInp(ppxIdxInp);
        ds        = sampleList.var4{itr};
        dsInp     = round(ds*srRatio);

        % Interpolate signals
        srawInp = spline(meta.t,S.raw,tInp);
        saccInp = spline(meta.t,S.acc,tInp);
        svelInp = spline(meta.t,S.vel,tInp);
        sdspInp = spline(meta.t,S.dsp,tInp);
        % clf; hold on; plot(meta.t,S.dsp,'-xr'), plot(tInp,sdspInp,'-dk');
        % set(gca,'xlim',[tppxOrg-.06 tppxOrg+.08])
        % plot(tppxOrg,0,'dr','lineWidth',2,'markerSize',11)
        % plot(tppxInp,0,'sk','lineWidth',2,'markerSize',11)
        
        % Equalise displacement polarities
        if params.equalPolarities
            diPolarity = params.diPolarity;   % No. of samples at which polarity is measured;
            
            if sdspInp(ppxIdxInp+diPolarity)<0;
            	srawInp = -srawInp;
                saccInp = -saccInp;
                svelInp = -svelInp;
                sdspInp = -sdspInp;
            end
        end
    
            
        % Save meta info
        META.m(itr)         = sampleList.m(itr);
        META.r(itr)         = sampleList.hypDist(itr);
        META.sr(itr)        = sampleList.sRate(itr);
        META.tppx(itr)      = tppxOrg;
        META.ppxIdx{itr}(1) = ppxIdxInp;
        META.ppxIdx{itr}(2) = dsInp;
        
        % Write vectors to cell-arrays
        TIME{itr} = tInp;
        RAW{itr}  = srawInp;
        ACC{itr}  = saccInp;
        VEL{itr}  = svelInp;
        DSP{itr}  = sdspInp;
        
        %         % Compare original and interpolated-truncated signals
        %         figure(1); clf; hold on;
        %         plot(tOrg   ,S.dsp,'or')
        %         plot(TIME{itr},DSP{itr},'-xk')
        %         plot(TIME{itr}(ppxIdxTrun),0,'sk','lineWidth',3)
        %         %set(gca,'xlim',[tppxOrg-.05, tppxOrg+.1])
        %         pause
    end
    
    
    % Save arrays to mat-file
    out.ACC  = ACC;
    out.VEL  = VEL;
    out.DSP  = DSP;
    out.TIME = TIME;
    out.META = META;
    
    % Interpolation information
    INP.dtInp = dtInp;
    INP.diUp  = diUp;
    INP.diLo  = diLo;
    out.INP   = INP;
    
    
    if params.saveOut
        
        outFileFullName = '/Users/mameier/programs/filterBank/var/wformMat/new/wform_onsets_interp100sps.mat';
        VAR.outFileFullName = outFileFullName;
        out.VAR             = VAR;
        
        fprintf(1,'Existing files will be overwritten... check if thats ok, or press enter to proceed.\n')
        pause
        save(outFileFullName,'out')
    end