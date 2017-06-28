function set_all_overwritten_fields_to_zero(tmpList)

names = fieldnames(tmpList);
nf    = size(names,1);
ntr   = numel(tmpList.(names{1}));
for i = 1:nf
    nentries = numel(tmpList.(names{i}));
    if nentries~=ntr
        tmpList.(names{i}) = 0;
        fprintf(1,sprintf('Note: Field %s only has %i entries; was set to a single 0\n',names{i},nentries))
    end
    
end
