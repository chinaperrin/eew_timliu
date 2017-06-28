function [tout,sout] = seismo2prelogram(t,s,tppx)

% Turns a seismogram into a prelogram as described by Kanamori and Mori 2000
sout = zeros(size(s));

% Truncate to positive time values
lgc  = logical(t>tppx);
sout = sout(lgc);
s    = s(lgc);
tout = t(lgc);
dt   = t(2)-t(1);
tout = tout-tout(1)+dt;

pos  = logical(s>0);
neg  = logical(s<0);
zilt = logical(s==0);

sout(pos)  = s(pos).^(1/3);
sout(neg)  = -(-s(neg)).^(1/3);
sout(zilt) = s(zilt);