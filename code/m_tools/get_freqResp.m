function [figHandle,Y,ff] = get_freqResp(A,B,sr,f0,figNum,o_newFig)

if (nargin<5)
    figNum   = 666;
    o_newFig = 1;    
end

if (o_newFig)
    col = 'k';
    lW  = 2;
else
    col = 'r';
    lW  = 2;
end



[Y ff]  = freqz(B,A,1024,sr);

figHandle = figure(figNum);
if (o_newFig==1); clf; end
hold on
plot(ff,abs(Y),'color',col,'lineWidth',lW)

grid on
set(gca,'xScale','log')
set(gca,'yScale','log')
set(gca,'xLim',[0.07,170])
[ylims] = get(gca,'ylim');
line([f0 f0],[ylims(1) ylims(2)],'Color','b')
%legend(['f0 = ',num2str(f0)])
title('Frequency Response', 'fontSize',15)
xlabel('Frequency', 'fontSize',15)
ylabel('Amplitude', 'fontSize',15)

