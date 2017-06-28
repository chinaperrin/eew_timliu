function [snr,noiseAmps,signalAmps,flg_issue] = get_snr(s,ppxIdx,sr,signalWindow,noiseWindow)
% Measure noise and signal amplitudes as well as signal-to-noise ratio (SNR)
%
% ppxIdx       = pick-index that separates noise and signal windows
% signalWindow = window length for measuring signal after pick [sec] 
% noiseWindow  = window length for measuring noise before pick [sec] 
%
% menandrin@gmail.com, 150103

eIdx = ppxIdx+signalWindow*sr-1;     % End   index of signal window
sIdx = ppxIdx-noiseWindow*sr;        % Start index of noise  window
    
nnw = noiseWindow *sr;
nsw = signalWindow*sr;

%  If signal is not long enough before or after pick, raise flag
flg_issue = false;

%% PRE-amps aka noise
if (ppxIdx>nnw)
    noiseAmps = s(sIdx:ppxIdx-1);
else
    noiseAmps = s(1:ppxIdx-1);
    fprintf(1,'8ung: Noise window truncated.\n')
    flg_issue = true;
end 

% Use the absolute value of the more extreme percentile 
%noise = max(abs(prctile(noiseAmps,[pctl, 100-pctl])));

%% POST-amps aka signal
if (ppxIdx+nsw<numel(s))
    signalAmps = s(ppxIdx:eIdx);
else
    signalAmps = s(ppxIdx:end);
    fprintf(1,'8ung: Signal window truncated.\n')
    flg_issue = true;
end

% Use the absolute value of the more extreme percentile 
%signal = max(abs(prctile(signalAmps,[pctl, 100-pctl])));

%% Power / SNR
signalPower = 1/signalWindow*sum(signalAmps.^2);
noisePower  = 1/noiseWindow *sum(noiseAmps .^2);
snr         = signalPower/noisePower;
