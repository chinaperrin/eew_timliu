function [snr,noiseAmps,signalAmps,flg_issue] = get_snr(s,ppxIdx,sr,signalWindow,gapWindow,noiseWindow)
% Measure noise and signal amplitudes as well as signal-to-noise ratio (SNR)
%
% ppxIdx       = pick-index that separates noise and signal windows
% signalWindow = window length for measuring signal after pick [sec] 
% gapWindow    = window length between end of noiseWindow and pick [sec]
% noiseWindow  = window length for measuring noise before pick [sec] 
%
% menandrin@gmail.com, 150103

% ...................................................................................
% INDICES:
%                   ||  noiseW   || gapW  ||    signalW      ||     
%  -----------------  -----------  -------  -----------------  ----------------  ...
%
%                    |            |        |
% start-indices:   snIdx        sgIdx    ssIdx (=ppxIdx)
%
%                                |        |                  |
% end-indices:                 enIdx    egIdx              esIdx
% ...................................................................................

    
nsn = noiseWindow *sr;
nsg = gapWindow   *sr;
nss = signalWindow*sr;

%  If signal is not long enough before or after pick, raise flag
flg_issue = false;


%% PRE-amps aka noise
snIdx = ppxIdx - nsg - nsn - 1;
if snIdx<1
    snIdx = 1;
    fprintf(1,'8ung: No room for gap between noise and signal windows.\n')
    flg_issue = true;
end
enIdx = snIdx + nsn -1;

if enIdx >= ppxIdx
    enIdx = ppxIdx-1;
    fprintf(1,'8ung: Noise window truncated.\n')
    flg_issue = true;
end
noiseAmps = s(snIdx:enIdx);


%% POST-amps aka signal
ssIdx = ppxIdx;
esIdx = ssIdx + nss -1;
if esIdx>length(s)
    esIdx = length(s);
    fprintf(1,'8ung: Signal window truncated.\n')
    flg_issue = true;
end
signalAmps = s(ssIdx:esIdx);



%% Power / SNR
signalPower = 1/signalWindow*sum(signalAmps.^2);
noisePower  = 1/noiseWindow *sum(noiseAmps .^2);
snr         = signalPower/noisePower;