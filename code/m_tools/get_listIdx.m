function idx = get_listIdx(traceFullName,trList)

[recordName] = get_recordName(traceFullName);
hit          = regexp(trList.fullName,recordName,'match');
idx          = find(cellfun(@(x) ~isempty(x), hit));

if (numel(idx)>1)
    fprintf(1,'more than one match\n')
elseif isempty(idx)
    fprintf(1,'no match\n')
end