function [amax_snp,idx_amax_snp] = filterBank_fulltrace_TFR(s,sr,ns,fc,npad)
% Returns maximuim amplitudes in all frequency bands specified by fc, from
% first to last snippet.

% INPUTS 
% s      = signal/waveform
% sr     = sampling rate
% pxIdx  = index of pick-time
% ns     = number of samples per snippet
% fc     = passband corner frequencies


realisationType = 8;    % Filter implementation type
o_verbose       = false;

nbands       = size(fc,1);
nsnp         = length(s)/ns;
idx_snp      = 1:ns:nsnp*ns;

sout         = cell(nbands,1);
amax_snp     = zeros(nbands,nsnp);
idx_amax_snp = zeros(nbands,nsnp);


for iband = 1:nbands
    
    % * * * * * * * * * * * * * * * * * * * * * *~ * * * * * * * * * * *
    fLow        = fc(iband,1);    % Lower corner frequency (for hipass)
    fUp         = fc(iband,2);    % Upper corner frequency (for lopass)
    [S,~,~] = butter_pass_tdomain_f(s,fLow,fUp,sr,realisationType,npad,o_verbose);
    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    
    snpMatrix  = reshape(S,ns,nsnp);
    [amax,idx] = max(abs(snpMatrix),[],1);
    idx_amax   = idx+idx_snp-1;
    
%     % Test: plot waveform along with max-values
%     figure(1); clf; hold on
%     plot(S)
%     plot(idx_amax,S(idx_amax),'xr')
%     plot(idx_amax,amax,'ok')
    
    amax_snp(iband,:)     = amax;
    idx_amax_snp(iband,:) = idx_amax;
end