function [cellList,clusterList] = get_merged_cellList(zList,dm,dr,nmin,o_debug,s)

rng(s); %randi([1 6],20,1)  

% Merge data cells into clusters of cell until there are at least nmin data
% points in each cluster. 
% Strategy: Randomly select cell with ndat<nmin. If cell has unmerged 
%           neighbour cells, add randomly selected neighbour cells until 
%           there are enough data points in cluster. If there are no 
%           unmerged neighbours, find nearest neighbour in r-direction that 
%           belongs to a cluster & merge current cell to that cluster.

% menandrin@gmail.com, 140919


% Use 'ii' for linear indices, il for line- and ic for column-indices
%
% cList = 'current list'  , list of linear cell indices in current cluster
% nList = 'neighbour list', list of linear cell indices of neighbouring
%         cells
% ndat  = a nm-1 by nr-1 matrix, containing #data points per cell
%
% use 'line' instead of 'row' to avoid confusion with r=distance
% m-values as lines   --> il
% r-values as columns --> ic

% Outstanding tasks / ideas:
% - Voronoi tesselations instead of fixed cells
% - exclude R>50 for M<3
% - instead of random walk, directly find nearest neighbours that are non-zero


o_checks = 0;


%% 1. Compute number of data points in equally spaced m- & r-bins

r = zList.hypDist;
m = zList.m;

rVect = 0:dr:100;      % r cell corners
nr    = numel(rVect);

mVect = 2:dm:8;        % m cell corners
nm    = numel(mVect);

ndat = zeros(nm-1,nr-1);
for im = 1:nm-1
    mlo = mVect(im);
    mup = mVect(im+1);
    
    for ir = 1:nr-1
        rlo = rVect(ir);
        rup = rVect(ir+1);
        
        ndat(im,ir) = sum( (m<mup) & (m>=mlo) & (r<rup) & (r>=rlo) );
    end
end

if o_debug

    ndatplot            = ndat;
    ndatplot(ndat>nmin) = nmin;
    figure(33); clf; hold on;
    imagesc(rVect(1:end-1)+dr/2,mVect(1:end-1)+dm/2,ndatplot)
    colorbar
    set(gca,'xlim',[min(rVect) max(rVect)],'ylim',[min(mVect) max(mVect)])
end
%plot(50,2,'xw','lineWidth',2)
ndat_orig = ndat;



%% 2. Merge cells until all clusters contain at least nmin data points 

% Store linear indices of cells and corresponding #data points <ndat>
iimax = numel(ndat);

% Each of the <nm-1> * <nr-1> cells has an entry in cellList
cellList.cid      = zeros(iimax,1,'uint16');                    % Cluster identification number
cellList.isMerged = false(iimax,1);                             % Flag for whether cell is part of cluster or not

% Each cluster has a list of linear cell indices in clusterList 
clusterList.ii    = cell (iimax,1);                             % Linear cell indices of all cells in a given cluster
clusterList.rLines = cell(iimax,1);
clusterList.mLines = cell(iimax,1);

ni  = numel( find(ndat<nmin & ndat>0) ); 
cid = 0;    


while ni>0  % While there are cells with ndat<nmin

    
    % 1. Start new cluster: randomly pick 1 cell with ndat<nmin
    % ---------------------------------------------------------
    idx      = find(ndat<nmin & ndat>0);
    ii       = idx(randi([1 numel(idx)],1,1));
    n        = ndat(ii);
    
    % Use index of chosen cell to start list of cell indices for current 
    % cluster, and find indices of adjacent cells
    ii_cList = ii;            
    ii_nList = find_and_check_neighbourCells(ii,mVect,rVect,cellList.isMerged);

    if o_debug
        [il,ic] = ind2sub(size(ndat),ii);
        plot(rVect(ic)+dr/2,mVect(il)+dm/2,'dw','markerFaceColor','r','markerSize',12,'markerEdgeColor','w','lineWidth',2)
        pause(.01)
    end
    
    
    % 2. Fill cluster until n>=nmin
    % ----------------------------- 
    while n<nmin
    
        
        % CASE 1: there are unmerged neighbours; add them one by one
        if ~isempty(ii_nList)
            
            % Add randomly chosen neighbour to list of cells in current cluster
            iix                         = ii_nList(randi([1 numel(ii_nList)],1,1));
            ii_cList                    = [ii_cList, iix];
            n                           = n+ndat(iix);                     
            ndat(ii_cList)              = n;
            cellList.isMerged(ii_cList) = true;
            
            % Remove added cell from neighbourList, find its own neighbours 
            % and add them to neighbourList
            ii_nList         = ii_nList(ii_nList~=iix);
            ii_newNeighbours = find_and_check_neighbourCells(iix,mVect,rVect,cellList.isMerged);
            ii_nList         = unique([ii_nList; ii_newNeighbours]);
            
            if o_debug
                % Update colormap
                ndatplot            = ndat;
                ndatplot(ndat>nmin) = nmin;
                gcf; imagesc(rVect(1:end-1)+dr/2,mVect(1:end-1)+dm/2,ndatplot)
                colorbar; set(gca,'xlim',[min(rVect) max(rVect)],'ylim',[min(mVect) max(mVect)])
                
                % Plot cells
                [iln,icn] = ind2sub(size(ndat),ii_nList);   % Line- & column-indices of neighbours
                [ilm,icm] = ind2sub(size(ndat),ii_cList);   % Line- & column-indices of cells in cluster
                plot(rVect(icn)+dr/2,mVect(iln)+dm/2,'ow')
                plot(rVect(icm)+dr/2,mVect(ilm)+dm/2,'ow','markerFaceColor','y')
                1+1;
            end
            
            isNewCluster = true;
            
            
            
        % CASE 2: there are no unmerged neighbours; add current cluster with existing cluster
        else 
            
            ii      = ii_cList(1);
            [il,ic] = ind2sub(size(ndat),ii);
            
            nnewcells = numel(ii_cList);
            
            itmp  = [(1:nr-1)',repmat(il,nr-1,1)];            % col and line indices on same line
            iitmp = sub2ind(size(ndat),itmp(:,2),itmp(:,1));  % corresponding linear indices
            
            % Find the nearest cell (with same magnitude) that has ndat>nmin
            hasData  = ndat(iitmp)>=nmin;
            ivect    = 1:numel(hasData);
            ihasData = ivect(hasData);
            [~,idx]  = min(abs(ihasData-ic));
            icn      = ihasData(idx);                         % Column index of nearest neighbour with enough data
            iin      = sub2ind(size(ndat),il,icn);            % Linear index of nearest neighbour with enough data
        
            cid_tmp  = cellList.cid(iin);                     % Find cluster id of that neighbour
            ii_cList = [clusterList.ii{cid_tmp}, ii_cList];   % Load linear cell indices & add current linear index list
            
            % Count #data points after adding current field
            n_tmp          = ndat(ii_cList(1));               % #data points before adding current cells
            n              = n + n_tmp;                       % #data points that has been computed for 
            ndat(ii_cList) = n;                               % current cluster either under 1. or under 2./CASE 1.

            cellList.isMerged(ii_cList) = true;
            isNewCluster                = false;
        end
    end
    
    
    % 3. Update cellList & clusterList with info from finished cluster
    % ----------------------------------------------------------------

    % Collect summary info
    nfields = numel(ii_cList);
    [il,ic] = ind2sub(size(ndat),ii_cList);
    mtmp    = mVect(il);
    rtmp    = rVect(ic);

    if isNewCluster
        
        % Update lists
        cid                    = cid + 1;
        clusterList.ii{cid}    = ii_cList;
        cellList.cid(ii_cList) = cid;
        
        % Print summary
        strng = sprintf('New cell: #%i, %i fields, m-range [%3.1f-%3.1f], r-range [%3.1fkm-%3.1fkm]\n', ...
                cid, nfields, min(mtmp), max(mtmp)+dm, min(rtmp), max(rtmp)+dr);
        fprintf(1,strng)
        1+1;
        
        % Compute list of bounding lines of cluster
        bdLines                 = get_cluster_boundary_lines;
        clusterList.rLines{cid} = bdLines.r;
        clusterList.mLines{cid} = bdLines.m;
        gcf; hold on
        if o_debug
            line(bdLines.r,bdLines.m,'color','w')
        end
    else
        
        % Update lists
        cellList.cid(ii_cList)  = cid_tmp;
        clusterList.ii{cid_tmp} = ii_cList;
        
        % Print summary
        strng = sprintf('          .. added %i fields to existing cluster #%i\n',nnewcells, cid_tmp);
        fprintf(1,strng)
        
        % Compute list of bounding lines of cluster
        bdLines                     = get_cluster_boundary_lines;
        clusterList.rLines{cid_tmp} = bdLines.r;
        clusterList.mLines{cid_tmp} = bdLines.m;
        gcf; hold on
        if o_debug
            line(bdLines.r,bdLines.m,'color','w')
        end
    end
    
    % Update ni
    ni = numel(find(ndat<nmin & ndat>0));
end
        
  

    
    
    
%% 3. All cells that needed merging into clusters have been merged. Now add single cells with n>=nmin
fprintf(1,'\nAdding single cells with enough data to cellList ...\n\n')

iiSingles = find(ndat(:)>=nmin & cellList.isMerged==0);
for i = 1:numel(iiSingles)
    ii                  = iiSingles(i);
    cid                 = cid+1;
    clusterList.ii{cid} = ii;
    cellList.cid(ii)    = cid;
end

% Remove empty entries in clusterList
lastIdx            = find(cellfun(@(x) ~isempty(x), clusterList.ii),1,'last');
clusterList.ii     = clusterList.ii    (1:lastIdx);
clusterList.mLines = clusterList.mLines(1:lastIdx);
clusterList.rLines = clusterList.rLines(1:lastIdx);

clusterList.ndat      = ndat;
clusterList.ndat_orig = ndat_orig;
clusterList.mVect     = mVect;
clusterList.rVect     = rVect;


% Consistency checks
if o_checks
    
    % Check if every cluster has ndat>nmin
    ncl = numel(clusterList.ii);
    for icl = 1:ncl
        ii   = clusterList.ii{icl};
        nsum = sum(ndat_orig(ii));
        nn   = unique(ndat(ii));
        if numel(nn)~=1 || ~isequal(nsum,nn)
            fprintf(1,'8ung: #data points cant be right. Check. Now.\n')
        end
    end
    
    % Is nCells+nNaughts equal to numel(ndat(:))?
    nCells      = numel(cell2mat(clusterList.ii'));
    nNaughts    = sum(ndat(:)==0 & cellList.isMerged==false);
    discrepancy = nCells+nNaughts - numel(ndat(:));
    if discrepancy~=0
        fprintf(1,'8ung: some numbers dont add up the way they should. Check. Now.\n')
    end
end









    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ii_newNeighbours] = find_and_check_neighbourCells(ii,linVect,colVect,isMerged)
        % ii are linear indices, as opposed to row and column indices
        
        nl    = numel(linVect);
        nc    = numel(colVect);
        iimax = (nl-1)*(nc-1);
        
        % Magnitude of ii^th cell
        [a,b] = ind2sub([nl-1 nc-1],ii);
        mii   = linVect(a);
        
        
        % Linear indices of 4 possible neighbours
        ii_newNeighbours = [ii-1; ii+1; ii-nl+1; ii+nl-1];
        
        % Avoid cells with indices out of box
        ii_newNeighbours = ii_newNeighbours( ii_newNeighbours<=iimax & ii_newNeighbours>0 );
        
        % Avoid cells that have already been merged with other cells
        ii_newNeighbours = ii_newNeighbours( ~isMerged(ii_newNeighbours));
        
        % Avoid comining cells from opposite side of parameter space
        [a,~]            = ind2sub([nl-1 nc-1],ii_newNeighbours);
        mnew             = linVect(a);
        mDiff            = abs(mii-mnew);
        ii_newNeighbours = ii_newNeighbours(mDiff<3);
        
    end


    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function bdLines = get_cluster_boundary_lines
    
        % 1. Make a list of 4 surrounding lines from each cell
        rList = [];
        mList = [];
        for j = 1:numel(ii_cList)
            
            [a,b] = ind2sub(size(ndat),ii_cList(j));
            mlo   = mVect(a);
            mup   = mlo+dm;
            rlo   = rVect(b);
            rup   = rlo+dr;
           
            rList = [rList; rlo rup; rlo rup; rlo rlo; rup rup];
            mList = [mList; mlo mlo; mup mup; mlo mup; mlo mup];
        end
        
        % 2. Remove all lines that appear more than once
        nl     = size(rList,1);
        rbList = zeros(size(rList));    % Initiate lists that contain only bounding lines
        mbList = zeros(size(mList));
        cl     = 0;                     % Line count
         
        rList = roundn(rList,-1);
        mList = roundn(mList,-1);
        for iln = 1:nl
            cLine    = [rList(iln,:), mList(iln,:)];
            rowFound = ismember([rList, mList],cLine,'rows');
            if sum(rowFound)==1
                cl = cl+1;
                rbList(cl,:) = cLine(1:2);
                mbList(cl,:) = cLine(3:4);
            end
        end
        ilast     = find(rbList(:,1)~=0,1,'last');
        rbList    = rbList(1:ilast,:);
        mbList    = mbList(1:ilast,:);
        bdLines.r = rbList';
        bdLines.m = mbList';
    end
end