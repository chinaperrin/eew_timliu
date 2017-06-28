function [blackList] = add2blackList(idx_tr,trList,comment,blackListName)

blackListFullName = strcat(['~/programs/filterBank/var/blackLists/',blackListName]);
global iN

if (~exist(blackListFullName,'file'))
	fprintf(1,'\nCreating new blackList\n')
    blackList = traceList(0);
else
    fprintf(1,'\nLoading existing blackList\n')
    load(blackListFullName);
end


% Check if trace is already contained in list
traceFullName = trList.fullName{idx_tr};
hits          = regexp(blackList.fullName,traceFullName);
idx_hits      = find(cellfun(@(x) ~isempty(x), hits));

if (~isempty(idx_hits))
    
    print_asciiart('asciiart/thumbsup.txt')
    fprintf(1,['\n\t.. trace is already contained in blackList, not added again\n\t      comment:',blackList.comment{idx_hits}])
    
else 

    fprintf(1,['\n\t.. adding ',traceFullName, '\n\t   to ',blackListFullName,'\n\t   comment: ',comment,'\n'])
    
    % Extract list entry from trList, add comment and add it to blackList
    tmpList            = trList.selectSubList(idx_tr);
    tmpList.comment{1} = sprintf('bl: %s',comment);
    tmpList.var4{1}    = sprintf('%s / i%i',date,iN);
    blackList.appendList(tmpList);
    
    % Save updated blackList
    save(blackListFullName,'blackList')
    fprintf(1,'\n\tupdated blackList has been saved.\n')
end