function [pj] = get_bivarPrior_GRsphere(b,zprime,o_plotPrior,o_printPrior)

% Compute bivariate joined & normalised prior, using Gutenberg Richter for 
% magnitude and the surface area of a sphere segement (i.e. assuming uniform
% eq-distribution in space) for distance.
%
% Input:    b       G-R b-value
%           zprime  depth of seismogenic zone
%                   & plotting options
%
% Output:   pj      bivariate log-likelihood of prior
%
% menandrin@gmail.com, 150103

global rr mm ftSize iN

[RR,MM] = meshgrid(rr,mm); 

prr = sphere_segment_surface_area(RR,zprime);   % Distance prior
pmm = 10.^(-b*MM);                              % Magnitude prior
pj  = prr.*pmm;                                 % Joined prior
pj  = pj./sum(pj(:));                           % Normalised joined prior
pj  = log(pj);


%% Plot
if o_plotPrior
    
    % Compute 1D functions
    pm = 10.^(-b*mm);
    pm = pm./sum(pm(:));
    pr = sphere_segment_surface_area(10.^rr,zprime);
    pr = pr./sum(pr(:));
    
    figure(222); clf;
    whitebg('k')
    ax1 = subplot(4,3,[1,2,4,5]);
    ax2 = subplot(4,3,[3,6]);
    ax3 = subplot(4,3,[10,11]);
    
    % Plot bivariate prior
    subplot(ax1)
    surf(rr,mm,pj); 
    view([0 90]);
    posbak = get(ax1,'position');
    hc = colorbar('location','southoutside');
    set(ax1,'position',posbak)
    title('Joined log-prior (G-R / sphere segment)','fontSize',ftSize)
    ylabel(['Magnitude'],'fontSize',ftSize)
    xlabel(['Hypocentral distance [km]'],'fontSize',ftSize)
    rlab = [5,10,20,50,100];
    set(gca,'xtick',log10(rlab),'xtickLabel',rlab,'fontSize',ftSize)

    cbPos = get(hc,'position');
    cbPos(2) = cbPos(2)-0.05;
    set(hc,'position',cbPos)
    xlabel(hc,'Prior log-likelihood')
    
    % Plot magnitude prior
    subplot(ax2);
    plot(mm,log(pm),'r','lineWidth',2)
    grid on;
    ylabel('L_{prior}(m)','fontSize',ftSize)
    view([-90 90])
    set(ax2,'XAxisLocation','top','fontSize',ftSize)
    xtl = get(ax2,'xTickLabel');
    xt  = get(ax2,'xTick');
    xt  = xt(1:end-1);
    xtl = xtl(1:end-1);
    ylm = get(ax2,'ylim');
    set(ax2,'XAxisLocation','top','fontSize',ftSize)
    xlabel('Magnitude','fontSize',ftSize)
    
    
    % Plot distance prior
    subplot(ax3)
    plot(rr,log(pr),'r','lineWidth',2)
    grid on;
    ylabel('L_{prior}(r)','fontSize',ftSize)
    set(gca,'xtick',log10(rlab),'xtickLabel',rlab,'fontSize',ftSize)
    xlabel(['Hypocentral Distance [km]'],'fontSize',ftSize)    

    pos    = get(ax3,'position');
    pos(2) = pos(2)+0.08;
    set(ax3,'position',pos)

    if o_printPrior
        figDir  = strcat(['~/programs/filterBank/fig/i',num2str(iN),'/prior/new/']);
        figName = strcat([figDir,'prior.eps']);
        
        set(gcf,'PaperPositionMode','auto')
        %set(gcf,'PaperPositionMode','manual')

        if exist(figName,'file'); unix(['rm ', figName]); end
        print('-depsc2',figName)
    end
end