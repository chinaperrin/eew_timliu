mm = [1, 2, 3, 4, 5, 6, 7, 8];
nm = numel(mm);
n  = 1e4;    % No. of random samples for computing durations

ftSize = 15;

% Rupture duration [sec] & length [km]
[D,L,c] = get_eqSourceDuration_WellsCoppersmith94(mm,n);


hf = figure(3); clf; 
for im = 1:nm
    
    % Plot LENGTH cdf
    h2      = subplot(2,1,2); hold on; box on; grid on
    [fl,xl] = ecdf(L(:,im));
    plot(xl,fl, 'lineWidth',2,'color','k')

    % Annotate
    txt = sprintf('M%i',mm(im));
    xt  = median(xl);
    yt  = .92;
    text(xt,yt,txt,'fontSize',ftSize,'BackgroundColor','w')
    
    
    
    % Plot DURATION cdf
    h1      = subplot(2,1,1); hold on; box on; grid on
    [fd,xd] = ecdf(D(:,im));
    plot(xd,fd, 'lineWidth',2,'color','k')
   
    % Annotate
    xt = median(xd);
    text(xt,yt,txt,'fontSize',ftSize,'BackgroundColor','w')

end
set(h1,'xscale','log', 'fontSize',ftSize)
set(h2,'xscale','log', 'fontSize',ftSize)

title(sprintf('Source scalings of Wells and Coppersmith, 1994, BSSA, using a=%4.2f & b=%4.2f, vr=unif[2.4 3.0]km/s ', ...
    c.a, c.b), 'fontSize',ftSize)
xlabel(h1,'Subsurface Rupture Length L [km]', 'fontSize',ftSize)
ylabel(h1,'p(L''<=L)', 'fontSize',ftSize)

xlabel(h2,'Rupture Duration D (=L/vr) [sec]', 'fontSize',ftSize)
ylabel(h2,'p(D''<=D)', 'fontSize',ftSize)


% set(h2,'xlim',[7e-3 1e3])
%set(gcf,'PaperPositionMode','auto')
%print('-dpng','~/programs/seismo/fig/i39/srcDuration/new/wellscoppersmith_duration')