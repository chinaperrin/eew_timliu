function [sout] = movingAverage(s,n)

ns   = numel(s);
sout = zeros(size(s));

for is = 1:ns
    
    interval = is-n:is+n;
    interval = interval(interval>=1 &interval<=ns);
    sout(is) = mean(s(interval));
end