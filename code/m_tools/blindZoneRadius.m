function [bzRadius] = blindZoneRadius(depth,time)

% Distance s-phase front has propagated at <time> seconds after origin time
%
% IN        depth       [m]
%           time        [sec]
% OUT       bzRadius    [m]

if depth<1e3
    fprintf('8ung: You specified depth<1000m. Make sure depth is given in meters.\n')
end

vs = 3400;  % [m/sec]

bzRadius = sqrt(time^2*vs^2 - depth^2);