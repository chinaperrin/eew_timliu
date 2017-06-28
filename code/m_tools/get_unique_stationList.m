function [stationList] = get_unique_stationList(trList)


dLat = 1e-3;    % ~100m
dLon = 1e-3;    % ~100m

ntr = numel(trList.eq.m);
stationList.lat     = 0;
stationList.lon     = 0;
stationList.vs30    = 0;
stationList.name    = cell(1,1);
stationList.name{1} = 'dummy';
% stationList.nw      = cell(1,1);
% stationList.nw  {1} = 'dummy';

for itr = 1:ntr
    print_iteration_numbers(itr,ntr,'tenthousands')
    
    stLon  = trList.station.lon(itr);
    stLat  = trList.station.lat(itr);
    stName = trList.station.name{itr};
    vs30   = trList.station.vs30(itr);
    %nwName = trList.networkName{itr};
    
    sameName   = logical(strcmp(stationList.name,stName));
    sameLat    = logical(abs(stationList.lat-stLat)<dLat);
    sameLon    = logical(abs(stationList.lon-stLon)<dLon);
    idxMatches = find(sameName&sameLat&sameLon);
    
    
    nmatches   = numel(idxMatches);
    if nmatches==0; 
    
        % Add station to list
        stationList.lat  = [stationList.lat ; stLat ];
        stationList.lon  = [stationList.lon ; stLon ];
        stationList.name = [stationList.name; stName];
        stationList.vs30 = [stationList.vs30; vs30  ];
        %stationList.nw   = [stationList.nw  ; nwName];
        
    elseif nmatches>1
        fprintf('More than one matching station found; picking first. check.\n')
    end
end

% Remove dummy entry (first entry)
stationList.lat  = stationList.lat(2:end);
stationList.lon  = stationList.lon(2:end);
stationList.name = stationList.name(2:end);
stationList.vs30 = stationList.vs30(2:end);