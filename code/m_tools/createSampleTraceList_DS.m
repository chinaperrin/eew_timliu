function sampleList = createSampleTraceList_DS(trList,nsamples)
%CREATE A TRACELIST WITH <nsamples> ENTRIES FROM EACH DATA SET

dsNameList = unique(trList.dataSetName);
nds        = numel(dsNameList);
sampleList = traceList(0);
for ids = 1:nds
    dsName    = dsNameList{ids};
    ns        = nsamples;
    idx       = find(strcmp(dsName,trList.dataSetName));
    nidx      = numel(idx);
    if nidx<ns; ns=nidx; end 
	idx       = datasample(idx,ns,'Replace',false);
    trList_tmp = trList.selectSubList(idx);
    sampleList.appendList(trList_tmp);
end