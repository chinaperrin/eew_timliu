clear all

tt = 0:.01:5;
p0 = 0;
Pexp = p0 + exp(tt);
Pquad = p0 + tt.^2;
P1p5  = p0 + tt.^1.5;
P1p3  = p0 + tt.^1.3;


figure(22); clf;

% lin scale
subplot(1,2,1); hold on; box on; grid on;
p1 = plot(tt,Pexp,'lineWidth',2);
p2 = plot(tt,Pquad,'k','lineWidth',2);
p3 = plot(tt,P1p3,'r','lineWidth',2);
p4 = plot(tt,P1p5,'m','lineWidth',2);

l1=legend('exp','quadratic','t^{1.3}','t^{1.5}');
set(l1,'fontSize',15,'location','northWest')



% log scale
subplot(1,2,2); hold on; box on; grid on;
plot(tt,Pexp,'lineWidth',2)
plot(tt,Pquad,'k','lineWidth',2)
plot(tt,P1p3,'r','lineWidth',2)
plot(tt,P1p5,'m','lineWidth',2)
set(gca,'yscale','log')

