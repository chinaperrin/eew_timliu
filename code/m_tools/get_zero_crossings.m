function [crossingIdx] = get_zero_crossings(s,o)

if nargin<2; o.plot=0; end

% Find zero crossings: when product between neighbouring samples is
% negative. <idx> is sample index before zero-crossing, i.e. zero-crossing
% occurs between idx and idx+1.
s1          = s(2:end);
s2          = s(1:end-1);
sprod       = s1.*s2;
crossingIdx = find(sprod<0);
% firstIdx   = idx(find(idx>ppxIdx,1,'first')); % Index of first zero-crossing after ppxIdx


if o.plot
    clf; hold on;
    plot(s,'r')
    %plot(idx,0,'dy')
    iMat = [crossingIdx, crossingIdx+1];
    zMat = zeros(size(iMat));
    line(iMat',zMat','color','c')
end