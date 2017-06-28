function [acc_int] = acc2vel(acc,sr,method)

% <method> = 1   --> integration with cumsum;
% <method> = 2   --> integration with cumtrapz;

if (nargin<3)
    method    = 2;
end


dt = 1/sr;

if (method==1)
    acc_int = cumsum(acc).*dt;
end

if (method==2)
    acc_int = cumtrapz(acc)*dt;
end

