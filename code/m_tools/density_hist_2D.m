function [weightMatrix,mEdges,dEdges] = get_weight_matrix(M,D,minimax,o_plot)

mmin = minimax.mmin;
mmax = minimax.mmax;
dmin = minimax.dmin;
dmax = minimax.dmax; 

if (nargin<4)
    o_plot = 1;
end

ftSize = 15;

% Magnitude
dm      = 0.1;                      % Sampling step
mEdges  = mmin-dm/2:dm:mmax+dm/2;   % Cell edges
mCenter = mEdges(1:end-1)+dm/2;     % Cell center
nm      = numel(mCenter);

% Distance
dd      = 1;
dEdges  = dmin-dd/2:dd:dmax+dd/2;
dCenter = dEdges(1:end-1)+dd/2;
nd      = numel(dCenter);

nVal = zeros(nm,nd);

for im = 1:nm
    
    % m-interval
    m_lo = mEdges(im);
    m_up = mEdges(im+1);
    
    for id = 1:nd
        
        % d-interval
        d_lo = dEdges(id);
        d_up = dEdges(id+1);
        
        % count values within intervals
        lgc         = (M>=m_lo) & (M<m_up) & (D>=d_lo) & (D<d_up);
        nVal(im,id) = sum(lgc);
    end
    1+1;
end



maxval=max(max(nVal));
title(num2str(maxval))
 
 
if o_plot
     
     figure(229); clf
     surf(dCenter,mCenter,nVal)
     view([0 90])
     colormap pink   % hot pink
     set(gca,'xlim',[dmin dmax],'ylim',[mmin mmax])
     set(gca,'fontSize',ftSize)
     xlabel('Distance [km]','fontSize',ftSize)
     ylabel('Magnitude','fontSize',ftSize)
     hcb = colorbar; set(hcb,'fontSize',ftSize)
     ylabel(hcb,'Counts','fontSize',ftSize)
     
 end
