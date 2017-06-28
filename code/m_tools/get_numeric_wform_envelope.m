function [env,idx] = get_numeric_wform_envelope(s,nsWidth,nsShift)

ns = numel(s);
iS = 1-nsShift;       % Make window start at one in first loop iteration

env    = zeros(ns,1);
idx    = zeros(ns,1);
icount = 0;
iE     = 0;
while iE<(ns-nsWidth-1)
    icount = icount + 1;
    
    iS          = iS+nsShift;
    iE          = iS+nsWidth-1;
    env(icount) = max(abs(s(iS:iE))); 
    idx(icount) = iS; 
end
env = env(1:icount);
idx = idx(1:icount)+(round(nsWidth/2));