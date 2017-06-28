function [y] = taper(y,dt,nt)

% Mulitplies initial nt samples of function y with a cosine(0:pi/2)
% so that strong oscillations during this initial period are suppressed,
% resulting in a smooth/gradual signal start.

% nt: number of samples over which taper is active
% dt: sampling interval

t1  = dt*[0:nt-1];
tap = sin(pi/2*1/t1(nt).*t1).^2;
tap = [tap ones(1,numel(y)-numel(tap))]';
y   = y.*tap;
