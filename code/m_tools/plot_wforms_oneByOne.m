function hf = plot_wforms_oneByOne(trList,fig)

    ntr = numel(trList.m);
    
    % PICKER (see SBPx.m for meaning of parameters)
    px            = load_picker_settings('production_short');
    fLow_prefilt  = px.Param.fLow_prefilt;
    fLow_px       = px.Param.fLow_px;
    fUp_px        = px.Param.fUp_px;
    ntap          = px.Param.ntap;
    if fig.plotPx; px.Opt.plotPx = 1; end
    
    % OPTIONS
    opt.process      = true;
    opt.fLow_prefilt = fLow_prefilt;
    opt.fOrder       = 2;
    opt.fMode        = 'causal';
    opt.ntap         = ntap;
    opt.intMode      = 'afterPx';
    opt.scp_wforms   = true;
    
    for itr = 1:ntr
        
        trList.printSingleTraceSummary(itr)
        fprintf(1,sprintf('%i/%i\tSNR: %e\n',itr,ntr,trList.snr(itr)))
        
        % VERTICAL: Read/process waveform, read meta-info
        [S,meta] = read_any_trace(trList.fullName{itr},trList,opt);
        t        = meta.t;
        ns       = numel(S.raw);
        sr       = meta.sr;
        instCode = trList.instrCode{itr};
        isSM     = (strcmp(instCode,'L') ||strcmp(instCode,'N') ||strcmp(instCode,'G'));
        ppxIdx   = trList.ppxIdx(itr);
        spxIdx   = trList.spxIdx(itr);
        tppx     = t(ppxIdx);
        if spxIdx<=ns; tspx = t(spxIdx);
        else           tspx = [];
        end
        pga      = trList.pga(itr);
        pgv      = trList.pgv(itr);
        pgd      = trList.pgd(itr);
        tpga     = t(trList.pgaIdx(itr)+ppxIdx);
        tpgv     = t(trList.pgvIdx(itr)+ppxIdx);
        tpgd     = t(trList.pgdIdx(itr)+ppxIdx);
        
        hf=figure(20); clf;
        %h = plot_wform(S.acc,t,tppx,tspx,[tppx-2 tppx+5],20);
        h = plot_wform(S.vel,t,tppx,tspx,[tppx-2 tppx+5],20);
        %h = plot_wform(S.dsp,t,tppx,tspx,[tppx-2 tppx+5],20);
        subplot(h.ax2); hold on;
        plot([tpga; tpga],[pga; -pga],'oy','lineWidth',2)
        plot([tpgv; tpgv],[pgv; -pgv],'sm','lineWidth',2)
        plot([tpgd; tpgd],[pgd; -pgd],'dc','lineWidth',2)
    end