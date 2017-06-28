function [sout] = filterBank_minimal(s,sr,fc,npad)
% Returns filtered waveforms from filterBank

% s  = signal/waveform
% sr = sampling rate
% fc = passband corner frequencies

realisationType = 8;    % Filter implementation type
o_verbose       = false;
nband           = size(fc,1);

sout = cell(nband,1);

for iband = 1:nband
    % * * * * * * * * * * * * * * * * * * * * * *~ * * * * * * * * * * *
    fLow        = fc(iband,1);    % Lower corner frequency (for hipass)
    fUp         = fc(iband,2);    % Upper corner frequency (for lopass)
    [sout{iband},~,~] = ...
        butter_pass_tdomain_f(s,fLow,fUp,sr,realisationType,npad,o_verbose);
    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
end
