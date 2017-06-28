function stationList = compile_unified_stationList(trList,optVs)


%% %%%%%%%%%%%
% NGA WEST 1 %
%%%%%%%%%%%%%%
fprintf('\nLooping over all NGA records to find unique list of stations ...')
idxNGA   = find(strcmp(trList.dataSetName,'ngawest1'));
ngawList = trList.selectSubList(idxNGA);
hasNgaw  = numel(ngawList.eq.m)~=0;
if hasNgaw
    stNgaw        = get_unique_stationList(ngawList);
    stNgaw.region = repmat({'ngawest1'},numel(stNgaw.name),1);
    %     fprintf('Getting VS30-values for all NGA West 1 sites ...')
    %     stNgaw.vs30   = get_vs30_ngawest1(stNgaw,ngawList,optVs);
else
    stNgaw.name   = [];
    stNgaw.lat    = [];
    stNgaw.lon    = [];
    stNgaw.region = [];
    stNgaw.vs30   = [];
end


%% %%%%%%
% SOCAL %
%%%%%%%%%
fprintf('\nLooping over all Californian records to find unique list of stations ...')
idxCali  = find(strcmp(trList.dataSetName,'scsn') |strcmp(trList.dataSetName,'scsnPx'));
caliList = trList.selectSubList(idxCali);
hasCali  = numel(caliList.eq.m)~=0;
if hasCali
    stCali   = get_unique_stationList(caliList);
    
    fprintf('Getting VS30-values for all Californian sites ...')
    [stCali.vs30,drVect] = get_vs30_california(stCali); 
    stCali.region = repmat({'cali'},numel(stCali.name),1);
else
    stCali.name   = [];
    stCali.lat    = [];
    stCali.lon    = [];
    stCali.region = [];
    stCali.vs30   = [];
end

%% %%%%%%
% JAPAN %
%%%%%%%%%
fprintf('\nLooping over all Japanese records to find unique list of stations ...')
idxJapan = find(strcmp(trList.dataSetName,'kNet') |strcmp(trList.dataSetName,'kikNet'));
japaList = trList.selectSubList(idxJapan);
hasJapa  = numel(japaList.eq.m)~=0; 
if hasJapa
    stJapan  = get_unique_stationList(japaList);
    fprintf('Getting VS30-values for all Japanese sites ...')
    stJapan.vs30 = get_vs30_japan(stJapan,optVs);
    stJapan.region = repmat({'japan'},numel(stJapan.name),1);
end



%% %%%%%%%
% OTHERS %
%%%%%%%%%%
isCali    = logical(strcmp(trList.dataSetName,'scsn') |strcmp(trList.dataSetName,'scsnPx'));
isJapa    = logical(strcmp(trList.dataSetName,'kNet') |strcmp(trList.dataSetName,'kikNet'));
isNGA     = logical(strcmp(trList.dataSetName,'ngawest1'));
idxOther  = find(~isCali &~isJapa &~isNGA);
otherList = trList.selectSubList(idxOther);
stOther   = get_unique_stationList(otherList);
stOther.vs30   = optVs.vs30_default*ones(size(stOther.lat,1),1); 
stOther.region = repmat({'other'},numel(stOther.name),1);




%% %%%%%%%%%%%%%
% FIND DOUBLES %
%%%%%%%%%%%%%%%%
% Find and eliminate stations that are in both Cali- and NGA-list
if hasNgaw &hasCali
    nst = numel(stNgaw.lat);
    hasCaliEntry = false(nst,1);
    for ist = 1:nst
        stLat  = stNgaw.lat (ist);
        stLon  = stNgaw.lon (ist);
        stName = stNgaw.name{ist};
        [idxStation,idxMultiples] = find_stationList_index(stLat,stLon,stName,stCali);
        
        if ~isempty(idxStation)
            hasCaliEntry(ist) = true;
        end
    end
    sum(hasCaliEntry)
    stNgaw.lat    = stNgaw.lat   (~hasCaliEntry);
    stNgaw.lon    = stNgaw.lon   (~hasCaliEntry);
    stNgaw.name   = stNgaw.name  (~hasCaliEntry);
    stNgaw.vs30   = stNgaw.vs30  (~hasCaliEntry);
    stNgaw.region = stNgaw.region(~hasCaliEntry);
end




%% %%%%%%
% UNIFY %
%%%%%%%%%
stationList.name   = [stCali.name  ; stJapan.name         ; stNgaw.name  ; stOther.name  ];
stationList.lat    = [stCali.lat   ; stJapan.lat          ; stNgaw.lat   ; stOther.lat   ];
stationList.lon    = [stCali.lon   ; stJapan.lon          ; stNgaw.lon   ; stOther.lon   ];
stationList.region = [stCali.region; stJapan.region       ; stNgaw.region; stOther.region];
stationList.vs30   = [stCali.vs30  ; stJapan.vs30.meanSlow; stNgaw.vs30  ; stOther.vs30  ]; 










% TTT .....................................................................
% nst = numel(stationList.nw);
% stationList.ttt = cell(nst,1);

% Grid for SoCal  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ngp   = 1001;                     % Number of gridpoints in either direction
% lat   = linspace(31.8,44,ngp)';
% scLAT = repmat(lat,1,ngp);
% lon   = linspace(-125,-113,ngp);
% scLON = repmat(lon,ngp,1);
% 
% 
% % Grid for Japan  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ngp   = 1001;
% lat   = linspace(30 ,45 ,ngp)';
% lon   = linspace(128,145,ngp);
% jpLAT = repmat(lat,1  ,ngp);
% jpLON = repmat(lon,ngp,1  );
% 
% ngp_zoom = 101;   % Number of grid points of focus area in x- and y-
%                         % directions (has to be uneven!)
% 
% % Add travel time tables
% % stationList       = get_ttt(stationList,'japan',jpLAT,jpLON);
% % stationList       = get_ttt(stationList,'socal',scLAT,scLON);
% stationList.jplat = jpLAT(:,1);
% stationList.jplon = jpLON(1,:);
% stationList.sclat = scLAT(:,1);
% stationList.sclon = scLON(1,:);