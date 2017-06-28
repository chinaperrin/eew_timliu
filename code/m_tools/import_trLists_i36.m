function [outList,fc,snpLength] = import_trLists_i36(dataSetNames,listType)

global iN
% listType      'traceList' or 'skipList'

nds = size(dataSetNames,1);
fprintf(1,['  Loading ',num2str(nds),' ',listType,'s ... '])

% Initiate global traceList
% if iN==37; outList = traceList_i37(0);
% else       outList = traceList(0);
% end
outList = traceList(0);

if     strcmp(listType,'traceList'); fName = 'trList.mat';
elseif strcmp(listType,'skipList');  fName = 'skipList.mat';
end
    

for ids = 1:nds

    fprintf(1,[num2str(ids),'   '])
    
    listName = strcat(dataSetNames{ids},fName);
    outName  = strcat(dataSetNames{ids},'fbOut.mat');
    tmp      = load(listName);
    load(outName)
    
    if strcmp(listType,'traceList')
        outList.appendList(tmp.TraceList)
    elseif strcmp(listType,'skipList')
        outList.appendList(tmp.SkipList)
    end

    %if ~exist('out.fc','var'); fc = 99999; out.fc = 99999; end
        
    % Check if the same parameters have been used for all data sets
    if (ids ~=  1)
        if (~isequal(fc,out.fc) || ~isequal(snpLength,out.snpLength) )
            fprintf(1,['\n\t8UNG: The data sets were produced with different parameters,', ...
            ' do not combine them!\n'])
            pause
        end
    end
    
    fc        = out.fc;
    snpLength = out.snpLength;
end
fprintf(1,' done.\n')
