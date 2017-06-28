function [azi1,azi2] = dsp2azi_la05(Z,E,N,interval,alpha)

Z  = Z(interval);
E  = E(interval);
N  = N(interval);
ns = numel(Z);

theta1 = zeros(ns,1);
R_ze   = zeros(ns,1);
R_zn   = zeros(ns,1);

for is = 2:ns;
    R_ze(is)   = alpha*R_ze(is-1) + Z(is)*E(is);
    R_zn(is)   = alpha*R_zn(is-1) + Z(is)*N(is);
    theta1(is) = rad2deg(pi + atan(R_ze(is)/R_zn(is)));
end
azi1 = mean(theta1(2:end));
azi2 = theta1(end); % Just tryin out ...