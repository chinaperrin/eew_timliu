% Invert linear system for relation between m_i, log(r_i) and narrow-band
% peak amplitudes A_ij
%
% A = d.c   --> c = dinv.A


nmax = numel(zList.m);  
ntr  = 1e2;             % Number of traces selected for inversion
idx  = randi(n,ntr,1);

spList = zList.selectSubList(idx);
m      = spList.m;
r      = log10(spList.hypDist);
d      = [ones(ntr,1), m, r];

% Extract peak amps: A = <n x 9> 
snp = 6;
A   = log10(cell2mat(cellfun(@(x) x(:,snp)', spList.amax,'UniformOutput',0)));

% Invert for coefficient matrix c (for non-square matrices, Matlab computes
% generalised inverse automatically)
c = d\A;

% Forward problem
%d = A/c;
d = A*pinv(c);

wHat = d(:,1);
mHat = d(:,2);
rHat = 10.^d(:,3);

% Compare predicted vs. catalog values
mRsd = m - mHat;
rRsd = 10.^r - rHat;
wRsd = 1 - wHat;


hist(mRsd)
median(mRsd)
std(mRsd)
hist(rRsd)


% Errors in magnitude an "ones" correlate:
plot(mRsd,wRsd,'xk')
plot(rRsd,wRsd,'xk')