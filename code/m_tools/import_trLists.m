function [outList] = import_trLists(dataSetNames,listType)

global iN
% listType      'traceList' or 'skipList'

fprintf(1,'Note: contents of ''prop''-field from first list are used\n')
nds = size(dataSetNames,1);
fprintf(1,['  Loading ',num2str(nds),' ',listType,'s ... '])

% Initiate global traceList
outList = traceList(0);

if strcmp(listType,'traceList');    fName = 'trList.mat';
elseif strcmp(listType,'skipList'); fName = 'skipList.mat';
end
    

for ids = 1:nds

    fprintf(1,[num2str(ids),'   '])
    
    listName = strcat(dataSetNames{ids},fName);
    if exist(listName,'file'); 
         
        tmp = load(listName);
        if strcmp(listType,'traceList'); 
            
            outList.appendList(tmp.TraceList)
            
            % Check if the same parameters have been used for all data sets
            fc            = outList.prop.fc;
            snpLength     = outList.prop.snpLength;
            fc_new        = tmp.TraceList.prop.fc;
            snpLength_new = tmp.TraceList.prop.snpLength;
            if (~isequal(fc_new,fc) || ~isequal(snpLength_new,snpLength) )
                fprintf(1,sprintf('\n\t8UNG: The data sets were produced with different parameters,\n\t      ==> do not combine them!\n'))
                pause
            end
        elseif strcmp(listType,'skipList');  
            
            outList.appendList(tmp.SkipList)
        end
        
    %else fprintf(1,sprintf('WARNING: list not found, not added to skipList: %s\n',listName))
    else fprintf(1,sprintf('(not found, not added) '))
    end
end
fprintf(1,' done.\n')
