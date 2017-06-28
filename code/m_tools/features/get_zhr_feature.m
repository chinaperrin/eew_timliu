function zhr = get_zhr_feature(sz,sen,sr,snpLength,nsnp)

% sz    vertical signal
% sen   horizontal signal

ns           = min([numel(sz), numel(sen)]);
nspersnippet = snpLength*sr;            % No. of samples per snippet
nsnpMax      = floor(ns/nspersnippet);
if nsnp>nsnpMax; nsnp=nsnpMax; end
iE           = (1:nsnp)'*nspersnippet;  % Indices of End of Waveform snippet

zhr = zeros(nsnp,1);
    for isnp = 1:nsnp;
 
        interval  = 1:iE(isnp); 
        zhr(isnp) = max(abs(sz(interval)))/max(abs(sen(interval)));
    end
