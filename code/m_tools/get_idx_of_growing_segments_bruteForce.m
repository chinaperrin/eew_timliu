function [idxGS,idxFS] = get_idx_of_growing_segments_bruteForce(s)

% Returns array of sample indices of growth episodes of signal <s>. If 
% s(is)>s(is<1), is is added to <idxGS>, otherwise it is added to <idxFS>.
% This is a brute force approach. Should think about how this could be done
% more efficiently.

ns         = numel(s);
wasGrowing = false;
iseg_g     = 0;     % Segment indices
iseg_f     = 0;
idxGS      = [];    % Array containing indices of growing (GE) and flat 
idxFS      = [];    % (FE) segmentes
s(1)       = 0;     % Force first sample comparison s(1)~=s(2) to be unequal, to initiate arrays

for is = 2:ns;
    %fprintf(1,sprintf('%i/%i\n',is,ns))
    isGrowing = logical(s(is)~=s(is-1));
    
    if isGrowing && ~wasGrowing
        % Start new growing segment
        iseg_g        = iseg_g+1;
        idxGS{iseg_g} = is;
        %fprintf(1,sprintf('Initiating %ith growth segement\n',iseg_g))
        
    elseif isGrowing && wasGrowing
        % Add is to growing segment
        idxGS{iseg_g} = [idxGS{iseg_g}; is];
        
    elseif ~isGrowing && wasGrowing
        % Start new flat segment
        iseg_f        = iseg_f+1;
        idxFS{iseg_f} = is;
        %fprintf(1,sprintf('Initiating %ith flat segement\n',iseg_f))
        
    elseif ~isGrowing && ~wasGrowing
        % Add is to flat segment
        idxFS{iseg_f} = [idxFS{iseg_f}; is];
    else
        fprintf(1,'Should never have arrived here. Flaw in logic. Find it. Now. Do it.\n')
    end
    
    % Back up previous state
    wasGrowing = isGrowing;
    1+1;
end