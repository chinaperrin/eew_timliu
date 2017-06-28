function h = cf_signals(s1,s2,figNum,xlms)

if isempty(figNum); h = gcf;
else                h = figure(figNum);
end

ftSize = 15;
whitebg('k')

% Plot signals
subplot(2,1,1); hold on; box on; grid on;
p1 = plot(s1,'w');
p2 = plot(s2,':r');
p3 = plot(s1-s2,'y');
l1 = legend('signal 1','signal 2','difference');
set(l1,'fontSize',ftSize)
set(gca,'fontSize',ftSize)
if ~isempty(xlms), set(gca,'xlim',xlms); end

% Plot difference
subplot(2,1,2); hold on; box on; grid on;
p3 = plot(s1-s2,'y');
l2=legend('diffenrence');
set(l2,'fontSize',ftSize)
set(gca,'fontSize',ftSize)
if ~isempty(xlms), set(gca,'xlim',xlms); end

% Plot maximum
maxdiff     = max(s1-s2);
maxdiff_pct = maxdiff/max(s1);
txt         = sprintf('max. diff.: %e, --> ~%5.1f percent',maxdiff,maxdiff_pct);
ylms        = get(gca,'ylim');
y           = ylms(1)+.15*(ylms(2)-ylms(1)); 
xlms        = get(gca,'xlim');
x           = xlms(1)+.5*(xlms(2)-xlms(1));
text(x,y,txt,'fontSize',ftSize);