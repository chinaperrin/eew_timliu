function maxstep = get_maxstep_feature(s,sr,snpLength,nsnp)

ns           = numel(s);
nspersnippet = snpLength*sr;            % No. of samples per snippet
nsnpMax      = floor(ns/nspersnippet);
if nsnp>nsnpMax; nsnp=nsnpMax; end
iE           = (1:nsnp)'*nspersnippet;  % Indices of End of Waveform snippet

maxstep.val = zeros(nsnp,1);
maxstep.idx = zeros(nsnp,1);
for isnp = 1:nsnp;

    interval          = 1:iE(isnp);
    [val,idx]         = max(diff(s(interval)*sr));
    maxstep.val(isnp) = val;
    maxstep.idx(isnp) = idx;
end
