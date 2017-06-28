function [ssrte] = clean_out_ssrte(ssrte)
% Goes through all fields at all levels of a ssrte-structure and replaces
% empty cell arays and matrices with a single empty entry. Note: could be
% coded more efficiently with a while loop: while hasNextLevel, add field
% to total field name; once deepest level has been reached, check for type
% of entry and do replacement.

%ssrte=evalList.var1{1};
%ssrte_bak = ssrte;
structName = 'ssrte';
o.verbose  = false;

% Level 1  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  
nmList1 = fieldnames(ssrte);
nf1     = numel(nmList1);
for if1 = 1:nf1
    fdName1     = nmList1{if1};
    fdFullName1 = sprintf('%s.%s',structName,fdName1);
    subStrcut1  = ssrte.(fdName1);
    hasNextLevel = isstruct(subStrcut1);
    
    % Level 2  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
    if hasNextLevel;  nmList2 = fieldnames(subStrcut1);
                      nf2     = numel(nmList2);
    else              nf2     = 0;
    end
    
    for if2 = 1:nf2
        fdName2      = nmList2{if2};
        fdFullName2  = sprintf('%s.%s',fdFullName1,fdName2);
        subStrcut2   = subStrcut1.(fdName2);
        hasNextLevel = isstruct(subStrcut2);
    
        
        % Level 3  .  .  .  .  .  .  .  .  .  .  .  .  .  .
        if hasNextLevel;  nmList3 = fieldnames(subStrcut2);
                          nf3     = numel(nmList3);
        else              nf3     = 0;
        end
        
        ctEmpty = 0;
        for if3 = 1:nf3
            fdName3      = nmList3{if3};
            fdFullName3  = sprintf('%s.%s',fdFullName2,fdName3);
            subStruct3   = subStrcut2.(fdName3);
            hasNextLevel = isstruct(subStruct3);
            if o.verbose; fprintf(1,sprintf('%s\n',fdFullName3)); end
            
            % Replace all empty matrices or cells with an empty cell: []
            if ~hasNextLevel
                subStruct3;
                if iscell(subStruct3)
                    isEmptyMat = cellfun(@(x) ~isempty(x),subStruct3);
                    if sum(isEmptyMat(:))==0
                        eval(sprintf('%s=[];', fdFullName3)); 
                        ctEmpty = ctEmpty+1;
                    else
                        [subStruct3_single] = double2single_cellArray(subStruct3);
                        eval(sprintf('%s=subStruct3_single;',fdFullName3)); 
                    end
                elseif ismatrix(subStruct3) 
                    if sum(subStruct3(:))==0; eval(sprintf('%s=[];',fdFullName3)); 
                                              ctEmpty = ctEmpty+1;
                    else                      eval(sprintf('%s=single(%s);',fdFullName3,fdFullName3)); 
                    end
                end
            else
                fprintf(1,'Has Level 4!\n'); pause; 
            end
        end
        % If no sub-fields are populated, clear entire field 
        %if ctEmpty==nf3 &&nf3~=0; eval(sprintf('%s=[];',fdFullName2)); end
    end
end
% toc
% whos ssrte_bak ssrte
% save('ssrte','ssrte')
% save('ssrte_bak','ssrte_bak')