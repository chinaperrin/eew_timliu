function [zsampleList,esampleList,nsampleList] = createSampleTraceList_DS_ZEN(zList,eList,nList,nsamples)
%CREATE A TRACELIST WITH <nsamples> ENTRIES FROM EACH DATA SET

dsNameList = unique(zList.dataSetName);
nds        = numel(dsNameList);
zsampleList = traceList(0);
esampleList = traceList(0);
nsampleList = traceList(0);
for ids = 1:nds
    dsName    = dsNameList{ids};
    ns        = nsamples;
    idx       = find(strcmp(dsName,zList.dataSetName));
    nidx      = numel(idx);
    if nidx<ns; ns=nidx; end 
	idx       = datasample(idx,ns,'Replace',false);
    zList_tmp = zList.selectSubList(idx);
    eList_tmp = eList.selectSubList(idx);
    nList_tmp = nList.selectSubList(idx);
    zsampleList.appendList(zList_tmp);
    esampleList.appendList(eList_tmp);
    nsampleList.appendList(nList_tmp);
end