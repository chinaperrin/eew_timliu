function [clusterList] = get_median_amps_in_cellClusters(zList,hList,o_debug)

% o_debug = 1;        % plot ndat as colormap with updates
rndSeed  = rng;
hundreds = linspace(1e2,1e5,1e3);

%[cellList,clusterList] = get_merged_cellList(zList,.1,5,10,o_debug,rndSeed);
[cellList,clusterList] = get_merged_cellList(zList,.1,5,10,0      ,randi([1 100],1,1) );

% Checks
%lgc=(cellList.cid==0);
%sum(clusterList.ndat(lgc));
%sum(clusterList.ndat_orig(lgc));

% Compute nbPGV percentiles in all clusters
ndat  = clusterList.ndat_orig;
mVect = clusterList.mVect;
rVect = clusterList.rVect;

ncl              = numel(clusterList.ii);
clusterList.pz05 = cell(ncl,nsnpmax);
clusterList.pz50 = cell(ncl,nsnpmax);
clusterList.pz95 = cell(ncl,nsnpmax);
clusterList.ph05 = cell(ncl,nsnpmax);
clusterList.ph50 = cell(ncl,nsnpmax);
clusterList.ph95 = cell(ncl,nsnpmax);

clusterList.m_median = zeros(ncl,1,'single');
clusterList.r_median = zeros(ncl,1,'single');

for icl = 1:ncl
    
    if icl==1 || ismember(icl,hundreds)
        promptString = sprintf('Computing amplitude percentiles for %ith of %i clusters ...\n',icl, ncl);
        fprintf(1,promptString)
    end
    
    % Find all zList indices in icl^th cluster
    iiList  = clusterList.ii{icl};  % Linear indices for ndat-matrix
    idxList = [];                   % traceList indices
    
    ncell = numel(iiList);
    
    for icell=1:ncell
        
        % Find cell boundary
        [il,ic] = ind2sub(size(ndat),iiList(icell));
        mlo     = mVect(il);
        mup     = mlo+(mVect(2)-mVect(1));
        rlo     = rVect(ic);
        rup     = rlo+(rVect(2)-rVect(1));
        
        % Find zList-indices of traces from this cell
        idx     = find( (m<mup) & (m>=mlo) & (r<rup) & (r>=rlo) );
        idxList = [idxList; idx];
    end
    %isequal(numel(idxList),unique(clusterList.ndat(iiList)))
    
    % Save indices to cluterList
    clusterList.idx{icl} = idxList;
    
    % Compute magnitude of cluster
    % - median?
    % - weighted mean?
    mList = zList.m(idxList);
    rList = zList.hypDist(idxList);
    
    clusterList.m_median(icl) = median(mList);
    clusterList.r_median(icl) = median(rList);
    
    % Compute vertical and horizontal median nbPGV
    for isnp = 1:nsnp
        
        snippet = snpList_ms(isnp);
        zAmps   = cell2mat(cellfun(@(x) x(:,snippet), zList.amax(idxList)','uniformOutput',0));
        hAmps   = cell2mat(cellfun(@(x) x(:,snippet), hList.amax(idxList)','uniformOutput',0));
        % clf; set(gca,'yscale','log'); grid on; hold on; plot(zAmps); plot(pz50,'lineWidth',4);plot(pz05,'lineWidth',2);plot(pz95,'lineWidth',2);
        
        clusterList.pz05{icl,snippet} = prctile(zAmps,5 ,2);
        clusterList.pz50{icl,snippet} = prctile(zAmps,50,2);
        clusterList.pz95{icl,snippet} = prctile(zAmps,95,2);
        
        clusterList.ph05{icl,snippet} = prctile(hAmps,5 ,2);
        clusterList.ph50{icl,snippet} = prctile(hAmps,50,2);
        clusterList.ph95{icl,snippet} = prctile(hAmps,95,2);
    end
end