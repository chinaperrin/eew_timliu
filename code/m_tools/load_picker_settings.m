function px = load_picker_settings(settingType)


clear px pxParam pxOpt pxWeight

if strcmp(settingType,'production_short') |isempty(settingType)
    pxParam.threshold    = 5;
    pxParam.searchRange  = [1 1];
    pxParam.nshift       = 1;
    pxParam.ntap         = 50;
    pxParam.fLow_prefilt = 0.075;
    pxParam.fLow_px      = 0.075;
    pxParam.fUp_px       = 30;
    pxParam.nsnpmin      = 1000;
    pxParam.signalWindow = .5;              % Window lenght for measuring signal amps & power after  p-pick [sec]
    pxParam.gapWindow    = .5;
    pxParam.noiseWindow  = .5;              % Window lenght for measuring noise  amps & power before p-pick [sec]
    pxParam.tw1_denom    = 2;               % tw1 = sr/tw1_denom
    pxParam.tw2_denom    = 10;              % tw2 = sr/tw2_denom
    pxWeight.wmax        = 1e1;
    pxWeight.wmin        = .1;
    pxWeight.signalWtFct = 'lin';
    pxWeight.noiseWtFct  = 'lin';
    pxOpt.plotPx         = 0;
    pxOpt.animatePx      = 0;
    pxOpt.quiet          = 1;

elseif strcmp(settingType,'production_long')
    pxParam.threshold    = 5;
    pxParam.searchRange  = [1 1];
    pxParam.nshift       = 1;
    pxParam.ntap         = 100;
    pxParam.fLow_prefilt = 0.075;
    pxParam.fLow_px      = 0.075;
    pxParam.fUp_px       = 30;
    pxParam.nsnpmin      = 1000;
    pxParam.signalWindow = 1;              % Window lenght for measuring signal amps & power after  p-pick [sec]
    pxParam.gapWindow    = .5;
    pxParam.noiseWindow  = 1;              % Window lenght for measuring noise  amps & power before p-pick [sec]
    pxParam.tw1_denom    = 2;               % tw1 = sr/tw1_denom
    pxParam.tw2_denom    = 10;              % tw2 = sr/tw2_denom
    pxWeight.wmax        = 1e1;
    pxWeight.wmin        = .1;
    pxWeight.signalWtFct = 'lin';
    pxWeight.noiseWtFct  = 'lin';
    pxOpt.plotPx         = 0;
    pxOpt.animatePx      = 0;
    pxOpt.quiet          = 1;

elseif strcmp(settingType,'safod') |isempty(settingType)
    pxParam.threshold    = 5;
    pxParam.searchRange  = [1 1];
    pxParam.nshift       = 1;
    pxParam.ntap         = 100;
    pxParam.fLow_prefilt = 0.075;
    pxParam.fLow_px      = 15;
    pxParam.fUp_px       = 100;
    pxParam.nsnpmin      = 1000;
    pxParam.signalWindow = .25;              % Window lenght for measuring signal amps & power after  p-pick [sec]
    pxParam.gapWindow    = .25;
    pxParam.noiseWindow  = .25;              % Window lenght for measuring noise  amps & power before p-pick [sec]
    pxParam.tw1_denom    = 2;               % tw1 = sr/tw1_denom
    pxParam.tw2_denom    = 10;              % tw2 = sr/tw2_denom
    pxWeight.wmax        = 1e1;
    pxWeight.wmin        = .1;
    pxWeight.signalWtFct = 'lin';
    pxWeight.noiseWtFct  = 'lin';
    pxOpt.plotPx         = 0;
    pxOpt.animatePx      = 0;
    pxOpt.quiet          = 1;

elseif strcmp(settingType,'narrow_band')
    pxParam.threshold    = 5;
    pxParam.searchRange  = [1 1];
    pxParam.nshift       = 1;
    pxParam.ntap         = 100;
    pxParam.fLow_prefilt = 0.075;
    pxParam.fLow_px      = 3;
    pxParam.fUp_px       = 10;
    pxParam.nsnpmin      = 1000;
    pxParam.signalWindow = .5;              % Window lenght for measuring signal amps & power after  p-pick [sec]
    pxParam.gapWindow    = .5;
    pxParam.noiseWindow  = .5;              % Window lenght for measuring noise  amps & power before p-pick [sec]
    pxParam.tw1_denom    = 2;               % tw1 = sr/tw1_denom
    pxParam.tw2_denom    = 10;              % tw2 = sr/tw2_denom
    pxWeight.wmax        = 1e1;
    pxWeight.wmin        = .1;
    pxWeight.signalWtFct = 'lin';
    pxWeight.noiseWtFct  = 'lin';
    pxOpt.plotPx         = 0;
    pxOpt.animatePx      = 0;
    pxOpt.quiet          = 1;

elseif strcmp(settingType,'debug')
    pxParam.threshold    = 5;
    pxParam.searchRange  = [1 1];
    pxParam.nshift       = 1;
    pxParam.ntap         = 50;
    pxParam.fLow_prefilt = 0.075;
    pxParam.fLow_px      = 0.075;
    pxParam.fUp_px       = 30;
    pxParam.nsnpmin      = 1000;
    pxParam.signalWindow = .5;              % Window lenght for measuring signal amps & power after  p-pick [sec]
    pxParam.gapWindow    = .5;
    pxParam.noiseWindow  = .5;              % Window lenght for measuring noise  amps & power before p-pick [sec]
    pxParam.tw1_denom    = 2;               % tw1 = sr/tw1_denom
    pxParam.tw2_denom    = 10;              % tw2 = sr/tw2_denom
    pxWeight.wmax        = 1e1;
    pxWeight.wmin        = .1;
    pxWeight.signalWtFct = 'lin';
    pxWeight.noiseWtFct  = 'lin';
    pxOpt.plotPx         = 1;
    pxOpt.animatePx      = 1;
    pxOpt.quiet          = 0;
end


px.Param  = pxParam;
px.Weight = pxWeight;
px.Opt    = pxOpt;
