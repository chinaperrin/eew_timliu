% Recompute regression coeffients for single event to check if everything is correct.

i = 1939;
eventId      = eqs.eventId(i);              % Choose random event
idxTarget    = eqs.traceId{i};              % Find resp. trace indices
idxTraining  = setdiff(1:nz,idxTarget);

zList.fullName(idxTarget)

r  = zList.hypDist(idxTraining);
m  = zList.m      (idxTraining);
wt = W.w3         (idxTraining);
n  = numel(r);
Yz = zeros(n,nbands);
Ye = zeros(n,nbands);
Yn = zeros(n,nbands);

for iband = 1:nbands
    Yz(:,iband) = cellfun(@(x) x(iband,snippet), zList.amax(idxTraining));
    Ye(:,iband) = cellfun(@(x) x(iband,snippet), eList.amax(idxTraining));
    Yn(:,iband) = cellfun(@(x) x(iband,snippet), nList.amax(idxTraining));
end
Yh = (Ye+Yn)./2;

clear Ye Yn

% 1. Coeffs from parfor-function
[C1] = get_regression_coeffs_parfor(Yz,Yh,r,m,wt,'CPB81mod2');


zTraining = zList.selectSubList(idxTraining);
eTraining = eList.selectSubList(idxTraining);
nTraining = nList.selectSubList(idxTraining);

% 2. Coeffs from regular function
o.plotRegFit=1;
[C2] = get_regression_coeffs(zTraining,eTraining,nTraining,'CPB81mod2',snippet,wt);

isequal(cell2mat(C1.zCoeff'),cell2mat(C2.zCoeff'))
isequal(cell2mat(C1.hCoeff'),cell2mat(C2.hCoeff'))

% 3. Coeffs from saved values
coeffFileFullName = strcat([out.outDir,'eqCoeff_CPB81mod2_w3_k',num2str(k),'_',num2str(nz),'traces.mat']);
load(coeffFileFullName)
C3 = Ca2{i};

isequal(cell2mat(C1.zCoeff'),cell2mat(C3.zCoeff'))
isequal(cell2mat(C1.hCoeff'),cell2mat(C3.hCoeff'))