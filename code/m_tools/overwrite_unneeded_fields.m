function overwrite_unneeded_fields(trList,fieldList)

% Because traceList.m is an object, fields cannot be removed. But you can
% still empty the space by replacing the field with a single integer.

fprintf(1,'Full size\t\t\t')
[~] = trList.printObjectSize;
nfd = numel(fieldList);
for ifd = 1:nfd
    
    fprintf(1,['Removing field ',fieldList{ifd},'\t\t'])

    eval(['trList.',fieldList{ifd},' = -1*ones(1,1,''uint8'');'])
    trList.printObjectSize;
end