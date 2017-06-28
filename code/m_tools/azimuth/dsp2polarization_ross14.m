function [rect,phi,p,s] = dsp2polarisation(Z,E,N,interval)

ni     = numel(interval);
covMat = 1./ni*[N(interval)*N(interval)', N(interval)*E(interval)', N(interval)*Z(interval)'; ...
                E(interval)*N(interval)', E(interval)*E(interval)', E(interval)*Z(interval)'; ...
                Z(interval)*N(interval)', Z(interval)*E(interval)', Z(interval)*Z(interval)'];

[eigenVect,eigenVal] = eig(covMat);


% Degree of polarization
lambda = sort([eigenVal(1,1),eigenVal(2,2),eigenVal(3,3)],'descend');
rect   = 1 - (lambda(2)+lambda(3))/(2*lambda(1));   % rectilinearity
phi    = acos(eigenVect(1,1)); % Incidence angle
p      = rect*cos(phi);
s      = rect*(1-cos(phi));
