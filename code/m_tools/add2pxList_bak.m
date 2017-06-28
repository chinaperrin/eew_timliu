function pxList = add2pxList(traceFullName,ppxIdx,snr,comment,pxListName)

pxListFullName = strcat(['~/programs/filterBank/var/pxLists/',pxListName]);

if (~exist(pxListFullName,'file')); 
	fprintf(1,'\nCreating new pxList\n')
    pxList.fullNames = 'dummy entry';
    pxList.ppxIdx    = 999999;
    pxList.snr       = 999999;
    pxList.comment   = 'dummy entry';

else
    fprintf(1,'\nLoading existing pxList\n')
    load(pxListFullName);
end

if isempty(comment); comment='nocomment'; end

% Check if trace is already contained in list
hits     = regexp(pxList.fullNames,traceFullName);
%idx_hits = find(cellfun(@(x) ~isempty(x), hits));

if ~isempty(hits)
    
    print_asciiart('asciiart/thumbsup.txt')
    fprintf(1,'\n\t.. trace is already contained in pxList, trace not added again.\n')
    fprintf(1,'\t   delete existing entry if you want to overwrite it.\n')
    
else 

    fprintf(1,['\n\t.. adding ',traceFullName, '\n\t   to ',pxListName,'\n\t   comment: ',comment,'\n'])
    
    % Add pick value and comments to pxList
    pxList.fullNames = {pxList.fullNames; traceFullName};
    pxList.ppxIdx    = [pxList.ppxIdx   ; ppxIdx       ];
    pxList.snr       = [pxList.snr      ; snr          ];
    pxList.comment   = {pxList.comment  ; comment      };
    
    % Save updated blackList
    save(pxListFullName,'pxList')
end