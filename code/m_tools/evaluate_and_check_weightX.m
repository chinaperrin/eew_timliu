function evaluate_and_check_weightX(zList,eList,nList,snippet)

% Check if inferred X-factor for relative  weight of dlogr and dm is about
% correct. If so dY1 and dY2 should be ~ equal.
nz = numel(zList.m);

% Evaluate X for dM = X*dR for kNN
dummyWeight = ones(nz,1);
A           = get_regression_coeffs(zList,eList,nList,'simpleLin',snippet,dummyWeight);

for i = 1:9
    a = A.hCoeff{i}(1);
    b = A.hCoeff{i}(2);
    d = A.hCoeff{i}(3);
    X(i) = -b/d;
end
X = X(4); % X =  -0.2903   -0.3540   -0.3987   -0.4494   -0.5040   -0.5636   -0.6079   -0.6123   -0.6600
X = 0.5;
%X = 0.2;

i = 4;
a = A.hCoeff{i}(1);
b = A.hCoeff{i}(2);
d = A.hCoeff{i}(3);
X = -b/d;

m0       = 5;
r0       = 100;
r0prime  = log10(r0);
r0prime2 = r0prime+X;
r02      = 10^r0prime2;
[~,~,Ypred]   = lsq_simpleLin([a,b,d],r0 ,m0  ,1);       % Reference amplitude
[~,~,Ypred_m] = lsq_simpleLin([a,b,d],r0 ,m0+1,1);       % Amp at m0+1 
[~,~,Ypred_r] = lsq_simpleLin([a,b,d],r02,m0  ,1);       % Amp at r02 

dY1 = Ypred-Ypred_r;
dY2 = Ypred-Ypred_m;
[dY1,dY2]