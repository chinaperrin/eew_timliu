n=2^4;
[bandwidth,density,X,Y]=kde2d([m,d],n,[mmin,dmin],[mmax, dmax]);
%[bandwidth,density,X,Y]=kde2d([m,d]);

%D  = density./max(max(density));

maxval       = max(max(density));
weightMatrix = maxnr./maxval;
weightMatrix2 = ceil(maxnr./maxval);


mm = linspace(mmin,mmax,n);
dd = linspace(dmin,dmax,n);


figure(112); clf

subplot(2,1,1)
imagesc(mm,dd,flipud(log10(D)))
colorbar

subplot(2,1,2)
plot(m,d,'x')
colorbar