function [ppxIdx2,snr2] = repick_SBPx(fullName,trList,tppx0,px)


% Find, load and pre-process trace
[S,meta]    = read_any_trace(fullName,trList,1);
sr          = meta.sr;
[~,ppxIdx0] = min(abs(meta.t-tppx0));
spx         = bworth(S.vel,sr,[px.fLow_px,px.fUp_px],'band',px.fOrder,'causal');  

% Pick it
[ppxIdx2,snr2,~,~] = SBPx(spx,1/sr,ppxIdx0,px.pxThreshold,px.sWindow,px.nWindow, ...
    px.wmax,px.wmin,px.sWtFct,px.nWtFct,px.ntap,px.plotPx,px.animatePx);

% Optionally add it to pick list