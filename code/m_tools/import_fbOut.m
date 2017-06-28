function [TraceList,SkipList,fc,snpLength] = import_fbOut(dataSetNames)

% Initiate global traceList
TraceList = traceList(0);
SkipList  = traceList(0);
nds       = size(dataSetNames,1);
fprintf(1,['  Of ',num2str(nds),' data sets, loading  '])

for ids = 1:nds

    fprintf(1,[num2str(ids),'   '])
    
    outName = strcat(dataSetNames{ids},'fbOut.mat');
    load(outName)

    TraceList.appendList(out.traceList)
    SkipList.appendList(out.skipList)
    
    
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