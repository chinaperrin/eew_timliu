function vs30Vect = get_vs30_ngawest1(stNgaw,trList,optVs)

error('ngaWest-traceLists already have populated vs30 fields from wformProc\n')
dLat = 1e-4;    % ~10m
dLon = 1e-4;    % ~10m

nst      = numel(stNgaw.lat);
vs30Vect = zeros(nst,1);
for ist = 1:nst

    % For each station in list, find all matching traceList entries
    stLon  = stNgaw.lon (ist);
    stLat  = stNgaw.lat (ist);
    stName = stNgaw.name{ist};
    
    sameName   = logical(strcmp(stName,trList.station.name));
    sameLat    = logical(abs(stLat-trList.station.lat)<dLat);
    sameLon    = logical(abs(stLon-trList.station.lon)<dLon);
    idxMatches = find(sameName&sameLat&sameLon);
    
    nmatches   = numel(idxMatches);
    if nmatches>=1
        tmp = unique(cell2mat(trList.var5(idxMatches)));
        if numel(tmp)>1; error('conflicting Vs30-values'); end

        vs30Vect(ist) = tmp;
    else
        fprintf(1,'No matching stations found in traceList. check.\n')
        vs30Vect(ist) = optVs.vs30_default;
    end
end