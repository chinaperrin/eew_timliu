function [A_short] = trim_amax(A,nsnpmin)

% Save lots of space by cutting <A> matrix when values no longer change.
% Find index <idxLastChange> when latest band reaches it's max value.
% Don't cut to less than <nsnpmin> entries, though.
% Note: time must be in second dimension.

[nbands,nsnp] = size(A);

if nbands>nsnp
   fprintf(1,'Thats too many bands and not enough snippets.... check. now.\n') 
end

isMax         = (A==repmat(A(:,end),1,size(A,2)));
idxLastChange = find(sum(isMax,1)==nbands,1,'first');


% Normal case: if last change was after nsnpmin snippets 
if (idxLastChange>=nsnpmin)
    A_short = A(:,1:idxLastChange);            % ... cut after last change
    
else
    % Pathological case 1: if last change was before nsnpmin snippets, and
    % record is long enough, i.e. longer than nsnpmin snippets ...
    if (nsnp>=nsnpmin)
        idxLastChange = nsnpmin;
        A_short = A(:,1:idxLastChange); % ... cut after last change
    
	% Pathological case 2: if last change was before nsnpmin snippets, but
	% record is shorter than nsnpmin snippets
    else
        
        lastCol = A(:,end);                  % ... replicate last column to
        nshort  = nsnpmin-nsnp;              %     reach nsnpmin columns
        A_short = [A, repmat(lastCol,1,nshort)];
    end
    
end