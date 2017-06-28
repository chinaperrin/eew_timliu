function [weightMat] = get_weight_matrix_counts(R,M,o_plot)

% Counts the data points in non-overlapping M/R-cells, and computes a
% weight proportional to 1/counts.

global minimax rEdges mEdges

fprintf(1,'\nComputing weight matrix based on data density ... ')

mmin = minimax.mmin;
mmax = minimax.mmax;
rmin = minimax.rmin;
rmax = minimax.rmax; 

if (nargin<3)
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


% Loop over all cells
counts = zeros(nr,nm);

for im = 1:nm
    
    % m-interval
    m_lo = mEdges(im);
    m_up = mEdges(im+1);
    
    for ir = 1:nr
        
        % d-interval
        r_lo = rEdges(ir);
        r_up = rEdges(ir+1);
        
        % Count values within intervals
        lgc           = (M>=m_lo) & (M<m_up) & (R>=r_lo) & (R<r_up);
        
        % Save values with the index of the lower cell edge
        counts(ir,im) = sum(lgc);
    end
end

maxCounts           = max(counts(:));
counts((counts==0)) = maxCounts;                  % Setting cells without observations to zero
weightMat           = round(maxCounts./counts);   % will give them a weight of 1

if (min(min(weightMat))<1)
    fprintf(1,'8UNG: weight = zero assigned to some cells. Check rounding.\n')
end

if o_plot
    
    figure(228); clf
    
    s1 = subplot(3,1,1);
    plot(R,M,'.k','markerSize',2)
    set(gca,'xlim',[rmin rmax],'ylim',[mmin mmax],'fontSize',ftSize)
    ylabel('Magnitude','fontSize',ftSize)
    %title(['Max. nr. of counts = ',num2str(maxCounts)],'fontSize',ftSize)
    set(s1,'xtickLabel',{''})
    %set(gca,'ytick',mmin:1:mmax)
    
    s2 = subplot(3,1,2);
    
    imagesc(rEdges(1:end-1)+dr/2,mEdges(1:end-1)+dm/2,counts')
    axis xy
    colormap pink   % hot pink
    set(gca,'xlim',[rmin rmax],'ylim',[mmin mmax],'fontSize',ftSize)
    ylabel('Magnitude','fontSize',ftSize)
    hcb = colorbar; set(hcb,'fontSize',ftSize)
    ylabel(hcb,'Counts','fontSize',ftSize)
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