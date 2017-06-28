function [idxStation,idxMultiples] = find_stationList_index(stLat,stLon,stName,stationList)

dLat = 1e-3;    % ~100m
dLon = 1e-3;    % ~100m

sameName   = logical(strcmp(stationList.name,stName));
sameLat    = logical(abs(stationList.lat-stLat)<dLat);
sameLon    = logical(abs(stationList.lon-stLon)<dLon);
idxMatches = find(sameName &sameLat &sameLon);

nmatches     = numel(idxMatches);
idxMultiples = [];
if     nmatches==1; idxStation = idxMatches;
elseif nmatches==0; idxStation = []; fprintf('Station not found in stationList\n')
else
    idxStation   = idxMatches(1); fprintf('More than one matching station found; picking first. check.\n')
    idxMultiples = idxMatches;
end

