function isdble = find_double_entries_in_cell(inCell)

n      = numel(inCell);
isdble = false(n,1);

for i = 1:n
    
    entry   = inCell{i};
    whosOut = whos('entry');
    
    if strcmp(whosOut.class,'double')
        isdble(i) = true;
    elseif strcmp(whosOut.class,'single')
        isdble(i) = false;
    else
        fprintf(1,'class not found, whats going on?\n')
    end
end