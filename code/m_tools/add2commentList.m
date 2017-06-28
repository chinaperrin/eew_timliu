function add2commentList(idx_tr,trList,comment,commentListName)

commentListFullName = strcat(['~/programs/filterBank/var/commentLists/',commentListName]);

if (~exist(commentListFullName,'file'))
    commentList = traceList(0);
else
    load(commentListFullName);
end


% Check if trace is already contained in list
traceFullName = trList.fullName{idx_tr};
hits          = regexp(commentList.fullName,traceFullName);
idx_hits      = find(cellfun(@(x) ~isempty(x), hits));

if (~isempty(idx_hits))
    
    print_asciiart('asciiart/thumbsup.txt')
    fprintf(1,'\nTrace is already contained in commentList, comment not added\n')
    
else 

    fprintf(1,['\nAdding ',traceFullName, ' and corecs \n\t to ',commentListFullName,'\n'])
    
    % Extract list entry from trList, add comment and add it to commentList
    idx                = get_listIdx(traceFullName,trList.fullName);
    tmpList            = trList.selectSubList(idx);
    tmpList.comment(:) = {comment};
    commentList.appendList(tmpList);
    
    
    % Save updated commentList
    save(commentListFullName,'commentList')
end
