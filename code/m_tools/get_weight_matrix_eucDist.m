function [weightMat] = get_weight_matrix_eucDist(R,M,fudge,o_plot)

% Computes data density with the function dataDensity.m which adds the
% Euclidian (?) distance from all data points wrt/ to the target point.
% Then computes weights that weigth*density = max(density) for all points.

global minimax rEdges mEdges
fprintf(1,'\nComputing weight matrix based on the squared euclidean distance of all data points ... ')

mmin = minimax.mmin;
mmax = minimax.mmax;
rmin = minimax.rmin;
rmax = minimax.rmax;

if (nargin<4)
    o_plot = 1;
end

ftSize = 15;

% -------------------------------------------------------------------------
% mE=mEdges; rE=rEdges;                         ir =      (rlo-rmin)/dr +1;
%                                               ir = floor( r -rmin)/dr)+1;
%
%         mmin                                            mmax
%
%         mE1                 mE2               mEn-1      mEn
%    rE1   |-------------------|--- ... // ... ---|---------|     rmin
%                              .
%          .  mE1<= m < mE2    .
%          .  rE1<= r < rE2    .
%          .                   .
%    rE2   -   .   .   .   .   .            
%          .
%    ...   //
%          .
%    rEn   -                                                      rmax
%    
% -------------------------------------------------------------------------


% Cell edges and intervals
nm = numel(mEdges)-1;
nr = numel(rEdges)-1;
dr = rEdges(2)-rEdges(1);
dm = mEdges(2)-mEdges(1);

%[density]  = dataDensity(R,M,nr,nm,[rmin;rmax;mmin;mmax],fudge);
[density]  = dataDensity(M,R,nm,nr,[mmin;mmax;rmin;rmax;],fudge);
maxdensity = max(density(:));
weightMat  = ceil(maxdensity./density);

if o_plot
    
    figure(229); clf
    
    s1 = subplot(3,1,1);
    plot(R,M,'.k','markerSize',2)
    set(gca,'xlim',[rmin rmax],'ylim',[mmin mmax],'fontSize',ftSize)
    ylabel('Magnitude','fontSize',ftSize)
    %title(['Max. nr. of counts = ',num2str(maxCounts)],'fontSize',ftSize)
    set(s1,'xtickLabel',{''})
    %set(gca,'ytick',mmin:1:mmax)
    
    s2 = subplot(3,1,2);
    
    imagesc(rEdges(1:end-1)+dr/2,mEdges(1:end-1)+dm/2,density')
    axis xy
    colormap pink   % hot pink
    set(gca,'xlim',[rmin rmax],'ylim',[mmin mmax],'fontSize',ftSize)
    ylabel('Magnitude','fontSize',ftSize)
    hcb = colorbar; set(hcb,'fontSize',ftSize)
    ylabel(hcb,'Density','fontSize',ftSize)
    set(s2,'xtickLabel',{''})
    box on
    
    s3 = subplot(3,1,3);
    imagesc(rEdges(1:end-1)+dr/2,mEdges(1:end-1)+dm/2,weightMat')
    axis xy
    set(gca,'xlim',[rmin rmax],'ylim',[mmin mmax],'fontSize',ftSize)
    xlabel('Hypocentral Distance [km]','fontSize',ftSize)
    ylabel('Magnitude','fontSize',ftSize)
    hcb = colorbar; set(hcb,'fontSize',ftSize)
    ylabel(hcb,'Weight','fontSize',ftSize)
    %set(gca,'ytick',mmin:1:mmax)
    
    % Positions
    dh = 1/38; % 4,10,1,10,1,10,2
    pos = get(s2,'pos');
    w   = pos(3);
    h   = 10*dh;
    xo = pos(1);
    yo1 = 26*dh;
    yo2 = 15*dh;
    yo3 =  4*dh;
    set(s1,'pos',[xo, yo1, w, h])
    set(s2,'pos',[xo, yo2, w, h])
    set(s3,'pos',[xo, yo3, w, h])
    
    %set(gcf,'PaperPositionMode','auto')
    %print('-depsc2','../../fig/data/weight/weightMat_i25.eps');
    %print('-dpng','-r300','../../fig/data/weight/weightMat_i25.png');
end

fprintf(1,'done.\n')