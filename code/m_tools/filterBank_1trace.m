function [amax,amaxIdx,sout] = filterBank_1trace(s,sr,pxIdx,snpLength,nsnp,fc,npad)
% Returns maximuim amplitudes in all frequency-band specified by fc between
% the pick-time <tpx> and <tpx + snippetLength>

% INPUTS 
% s             = signal/waveform
% sr            = sampling rate
% pxIdx         = index of pick-time
% snippetLength = length of time window of interest in seconds
% fc            = passband corner frequencies
% fc_pre        = corner frequency for pre-filtering (high-pass)

realisationType = 8;    % Filter implementation type
o_verbose       = false;

nband   = size(fc,1);

sout    = cell(nband,1);
amax    = zeros(nband,nsnp);
amaxIdx = zeros(nband,nsnp);

for iband = 1:nband
    
    % * * * * * * * * * * * * * * * * * * * * * *~ * * * * * * * * * * *
    fLow        = fc(iband,1);    % Lower corner frequency (for hipass)
    fUp         = fc(iband,2);    % Upper corner frequency (for lopass)
    [sout{iband},~,~] = ...
        butter_pass_tdomain_f(s,fLow,fUp,sr,realisationType,npad,o_verbose);
    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    % Compute maxima in each of the <nInc> time snippets
    % Find time-indicex of the pick
    snpInt  = round(sr*snpLength); 
    lastIdx = pxIdx+snpInt:snpInt:pxIdx+nsnp*snpInt;
    
    for isnp = 1:nsnp
        S                              = sout{iband}(pxIdx:lastIdx(isnp));
        [amax(iband,isnp),amaxIdx_tmp] = max(abs(S));
        amaxIdx(iband,isnp)            = amaxIdx_tmp + pxIdx-1;
    end
end
