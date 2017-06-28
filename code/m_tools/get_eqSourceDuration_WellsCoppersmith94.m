function [D,RLD,coeffs] = get_eqSourceDuration_WellsCoppersmith94(M,n)

% Compute source duration distributions using the scaling D~L/vr and
% length/magnitude relationships of Wells and Coppersmith, 1994, BSSA.
% Randomly samples coefficients a & b, as well as rupture velocity (from a
% uniform distribution). Uses the coefficients from all rupture mechanisms. 
% Returns lengths in [km] and durations in [sec]
%
% menandrin@gmail.com, 150423

% n = 1e5;
% M = (2:.1:9)'
nm = numel(M);

if size(M,1)<size(M,2); M = M'; end     % Turn M into a <1 x n> vector

% Randomly sample coefficients a & b
a  = -2.44;
sa = 0.11;
b  = 0.59;
sb = 0.02;
ar = a + sa*randn(n,1);
br = b + sb*randn(n,1);

coeffs.a  = a;
coeffs.b  = b;
coeffs.sa = sa;
coeffs.sb = sb;

% Randomly sample rupture velocity from uniform distribution in interval 
% [vrLo vrUp]
vrUp = 3.0;
vrLo = 2.4;
vr   = vrLo + (vrUp-vrLo)*rand(n,1); 

% Replicate values for matrix computations
MM = repmat(M',n,1);
Ar = repmat(ar,1,nm);
Br = repmat(br,1,nm);
Vr = repmat(vr,1,nm);


RLD = 10.^(Ar + Br.*MM);        % Subsurface rupture length in [km]
D   = RLD./Vr;                  % Source duration
