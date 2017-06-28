function cellOut = harmonise_cell_lengths(cellIn,modeName)

% modeName      'zeropad':  Add zeros so that all lines are as long as
%                           longest one
%               'truncate': Truncate all lines to be as short as shortest
%                           one
%               'replast':  Replicate last value until all lines are as
%                           long as longest one
%
% menandrin@gmail.com, 150515

nl      = numel(cellIn);
cellOut = cell(size(cellIn));

if strcmp(modeName,'zeropad')
    1+1;
    n    = cellfun(@(x) numel(x), cellIn);
    nmax = max(n);

    for il = 1:nl
        cellOut{il} = [cellIn{il};zeros(nmax-n(il),1)];
    end    
    
elseif strcmp(modeName,'truncate')
    
elseif strcmp(modeName,'replast')
 
    %     nmax = floor(tRange(2)/snpLength);
    %     ntr  = numel(idx);
    %     S    = zeros(ntr,nmax);
    %     
    %     for itr = 1:ntr
    %         ns = length(amps{itr});
    %         if ns>=nmax; S(itr,:) = amps{itr}(1:nmax);                                  % (1)
    %         else         S(itr,:) = [amps{itr}, repmat(amps{itr}(end),1,nmax-ns)];      % (2)
    %         end
    %     end
end


