function [snr,noise] = get_noise(s,ppxIdx,sr,signalWindow,noiseWindow)
% Evaluate signal amplitudes in the 3sec (if available) before and after p-pick

% signalWindow = window length for measuring signal after pick [sec] 
% noiseWindow  = window length for measuring noise before pick [sec] 
dt = 1/sr;

eIdx = ppxIdx+signalWindow*sr-1;     % End   index of signal window
sIdx = ppxIdx-noiseWindow*sr;        % Start index of noise  window

    
nnw = noiseWindow*sr;
nsw = signalWindow*sr;

% PRE-amps aka noise
if (ppxIdx>nnw)
    noise = max(abs(s(sIdx:ppxIdx-1)))-abs(mean(s(sIdx:ppxIdx-1)));
else
    noise = max(abs(s(1:ppxIdx-1)))   -abs(mean(s(1:ppxIdx-1)));
end

% POST-amps
if (ppxIdx+nsw<numel(s))
    post_amps = max(abs(s(ppxIdx:eIdx)))-abs(mean(s(ppxIdx:eIdx)));
else
    post_amps = max(abs(s(ppxIdx:end))) -abs(mean(s(ppxIdx:end)));
end

% SNR
snr = post_amps/noise;