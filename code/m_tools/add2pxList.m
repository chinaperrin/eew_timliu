function [pxList] = add2pxList(idx_tr,corrList,pxListName)

% corrList is the updated trace list that includes the corrected picks,
% recomputed peak amps, etc.


pxListFullName = strcat(['~/programs/filterBank/var/pxLists/',pxListName]);

if (~exist(pxListFullName,'file'))
	fprintf(1,'\nCreating new pxList\n')
    pxList = traceList(0);
else
    fprintf(1,'\nLoading existing pxList\n')
    load(pxListFullName);
end


% Check if trace is already contained in list
traceFullName = corrList.fullName{idx_tr};
hits          = regexp(pxList.fullName,traceFullName);
idx_hits      = find(cellfun(@(x) ~isempty(x), hits));

if (~isempty(idx_hits))
    
    print_asciiart('asciiart/thumbsup.txt')
    fprintf(1,['\n\t.. trace is already contained in pxList, not added again\n\t      comment:',pxList.comment{idx_hits},'\n'])
    
else 

    fprintf(1,['\n\t.. adding ',traceFullName, '\n\t   to ',pxListFullName,'\n\t   comment: ',pxList.comment{idx_hits},'\n'])
    
    % Extract list entry from corrList, add comment and add it to pxList
    tmpList = corrList.selectSubList(idx_tr);
    pxList.appendList(tmpList);
    
    % Save updated pxList
    save(pxListFullName,'pxList')
    fprintf(1,'Updated pxList has been saved.\n')
end