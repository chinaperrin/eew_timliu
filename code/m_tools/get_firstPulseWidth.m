function [pulseWidth,firstIdx] = get_firstPulseWidth(s,ppxIdx,sr,o_plot)

% Measure first zero crossing after p-pick on velocity trace 
% --> corresponds to first maximum in displacement.

if nargin<4; o_plot=false; end

% Find zero crossings: when product between neighbouring samples is
% negative. <idx> is sample index before zero-crossing, i.e. zero-crossing
% occurs between idx and idx+1.
s1         = s(2:end);
s2         = s(1:end-1);
sprod      = s1.*s2;
idx        = find(sprod<0);
firstIdx   = idx(find(idx>ppxIdx,1,'first')); % Index of first zero-crossing after ppxIdx
pulseWidth = (firstIdx-ppxIdx)/sr;

if o_plot
    clf; hold on; grid on;
    plot(s)
    plot(idx,0,'xr')
    plot(firstIdx,0,'ok','markerFaceColor','r','markerSize',12)
    plot(ppxIdx  ,0,'dk','markerFaceColor','b','markerSize',12)
end