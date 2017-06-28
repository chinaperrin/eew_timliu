function cf_px_after_wformProc_SMBB(catList)

% Use colocated BB and SM records to check the effect that two integration
% methods have on the Pd(t) curves.

opt.process = 1;
opt.intMode = 'afterPx';

fMode = 'causal';
fOrder = 4;

fLow_prefilt = 0.075;
fLo          = 1;
fUp1         = 10;
fUp2         = 30;


ntap         = 100;


pxThreshold = 5;
fOrder      = 4;
sWindow     = 1;              % Window lenght for measuring signal amps after p-pick [sec]
nWindow     = 1;              % Window lenght for measuring noise before p-pick [sec]
wmax        = 1e1;
wmin        = .1;
sWtFct      = 'exp';          % 'lin' or 'exp' or '1-exp'
nWtFct      = 'exp';          % 'lin' or 'exp' or '1-exp'
plotPx      = false;
animatePx   = false;
    
fprintf(1,'Finding double records  ...')
matches = find_colocated_SM_BB_records(catList);

% Remove multiple entries
lgc1          = (cellfun(@(x) numel(x), matches.bbIdx)>1);
lgc2          = (cellfun(@(x) numel(x), matches.smIdx)>1);
lgc           = logical(lgc1 |lgc2);
matches.bbIdx = matches.bbIdx(~lgc);
matches.smIdx = matches.smIdx(~lgc);

% Throw out all records/matches with hypDist>25km
lgc           = catList.hypDist(cell2mat(matches.bbIdx))<=25;
matches.bbIdx = matches.bbIdx(lgc);
matches.smIdx = matches.smIdx(lgc);


%% Compare integration methods
nmatch = numel(matches.bbIdx);
for imatch = 1:nmatch
    
    fprintf(1,sprintf('%i/%i\n',imatch,nmatch))
    
    ibb = matches.bbIdx{imatch};
    ism = matches.smIdx{imatch};
    
    catList.printSingleTraceSummary(ibb)
    catList.printSingleTraceSummary(ism)
    
    
    % BROAD BAND  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
    [v,meta_v] = read_any_trace(catList.fullName{ibb},catList,opt);
    tv         = meta_v.t;
	ppxIdx_v   = meta_v.ppxIdx;
    tppx_v     = tv(ppxIdx_v);
    sr_v       = meta_v.sr;        
    
    vm          = v.raw - mean(v.raw(1:ppxIdx_v));
    vtap        = taper(vm,sr_v,ntap);
    vh          = bworth(vtap,sr_v,fLow_prefilt,'high',fOrder,fMode);
    vb1         = bworth(vtap,sr_v,[fLo fUp1],'band',fOrder,'causal');
    vb2         = bworth(vtap,sr_v,[fLo fUp2],'band',fOrder,'causal');

    
    % STRONG MOTION   . . . . . . . . . . . . . . . . . . . . . . . . . . .
    [a,meta_a] = read_any_trace(catList.fullName{ism},catList,opt);
    ta         = meta_a.t;
	ppxIdx_a   = meta_a.ppxIdx;
    tppx_a     = ta(ppxIdx_a);
    sr_a       = meta_a.sr;        
    
    am          = a.raw - mean(a.raw(1:ppxIdx_a));
    atap        = taper(am,sr_a,ntap);
    ab1         = bworth(atap,sr_a,[fLo fUp1],'band',fOrder,'causal');
    ab2         = bworth(atap,sr_a,[fLo fUp2],'band',fOrder,'causal');

    ah          = bworth(atap,sr_a,fLow_prefilt,'high',fOrder,fMode);
    ahi         = integrate_wform(ah,sr_a,ppxIdx_a,opt);                      % Integrate to vel (from [m/s^2] to [m/s])
    ahih        = bworth(ahi ,sr_a,fLow_prefilt,'high',fOrder,fMode);         % Highpass

    ahib1       = bworth(ahi,sr_a,[fLo fUp1],'band',fOrder,'causal');
    ahib2       = bworth(ahi,sr_a,[fLo fUp2],'band',fOrder,'causal');

    
    % PICK ALL VERSIONS   . . . . . . . . . . . . . . . . . . . . . . . . . 
    [pxv.vh   ,~,~,~] = SBPx(vh   ,1/sr_v,ppxIdx_v-150,pxThreshold,sWindow,nWindow,wmax,wmin,sWtFct,nWtFct,ntap,plotPx,animatePx);
    [pxv.vb1  ,~,~,~] = SBPx(vb1  ,1/sr_v,ppxIdx_v-150,pxThreshold,sWindow,nWindow,wmax,wmin,sWtFct,nWtFct,ntap,plotPx,animatePx);
    [pxv.vb2  ,~,~,~] = SBPx(vb2  ,1/sr_v,ppxIdx_v-150,pxThreshold,sWindow,nWindow,wmax,wmin,sWtFct,nWtFct,ntap,plotPx,animatePx);
    [pxa.ah   ,~,~,~] = SBPx(ah   ,1/sr_a,ppxIdx_a-150,pxThreshold,sWindow,nWindow,wmax,wmin,sWtFct,nWtFct,ntap,plotPx,animatePx);
    [pxa.ab1  ,~,~,~] = SBPx(ab1  ,1/sr_a,ppxIdx_a-150,pxThreshold,sWindow,nWindow,wmax,wmin,sWtFct,nWtFct,ntap,plotPx,animatePx);
    [pxa.ab2  ,~,~,~] = SBPx(ab2  ,1/sr_a,ppxIdx_a-150,pxThreshold,sWindow,nWindow,wmax,wmin,sWtFct,nWtFct,ntap,plotPx,animatePx);
    [pxa.ahih ,~,~,~] = SBPx(ahih ,1/sr_a,ppxIdx_a-150,pxThreshold,sWindow,nWindow,wmax,wmin,sWtFct,nWtFct,ntap,plotPx,animatePx);
    [pxa.ahib1,~,~,~] = SBPx(ahib1,1/sr_a,ppxIdx_a-150,pxThreshold,sWindow,nWindow,wmax,wmin,sWtFct,nWtFct,ntap,plotPx,animatePx);
    [pxa.ahib2,~,~,~] = SBPx(ahib2,1/sr_a,ppxIdx_a-150,pxThreshold,sWindow,nWindow,wmax,wmin,sWtFct,nWtFct,ntap,plotPx,animatePx);
    
    % Turn pick indices into pick times
    tppxv.vh    = tv(pxv.vh);
    tppxv.vb1   = tv(pxv.vb1);
    tppxv.vb2   = tv(pxv.vb2);
    
    tppxa.ah    = ta(pxa.ah);
    tppxa.ab1   = ta(pxa.ab1);
    tppxa.ab2   = ta(pxa.ab2);
    tppxa.ahih  = ta(pxa.ahih);
    tppxa.ahib1 = ta(pxa.ahib1);
    tppxa.ahib2 = ta(pxa.ahib2);
    
    
    sv.vh    = vh;
    sv.vb1   = vb1;
    sv.vb2   = vb2;
    sv.sraw  = v.raw;
    
    sa.ah    = ah;
    sa.ab1   = ab1;
    sa.ab2   = ab2;
    sa.ahih  = ahih;
    sa.ahib1 = ahib1;
    sa.ahib2 = ahib2;
    sa.sraw  = a.raw;
    
    plot_wform_components(sv,tv,tppxv,[tppx_v-.4 tppx_v+.8],18)
    plot_wform_components(sa,ta,tppxa,[tppx_a-.4 tppx_a+.8],19)
    1+1;
end