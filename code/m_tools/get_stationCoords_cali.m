function [stLat,stLon,nw,sta,chan] = get_stationCoords_cali(traceFullName,stationList)


flg_found = false;

ptIdx    = regexp(traceFullName,'\.');

nw   = traceFullName(ptIdx(1)+1:ptIdx(2)-1);
sta  = traceFullName(ptIdx(2)+1:ptIdx(3)-1);
chan = traceFullName(ptIdx(3)+1:ptIdx(4)-1);

% find station
lgc_nw  = strcmp(stationList.nw,nw);
lgc_sta = strcmp(stationList.name,sta);
idx     = find(lgc_nw.*lgc_sta);

if (~isempty(idx))
    
    % check if corresponding channel is in list
    match     = regexp(stationList.chan(idx),chan(1:2));
    idx_match = find(cellfun(@(x) ~isempty(x), match));
    
    if (~isempty(idx_match))
        
        idx       = idx(idx_match(1));    % If multiple matches, take first one
        flg_found = true;
        
        stLat = stationList.lat(idx);
        stLon = stationList.lon(idx);
        nw    = stationList.nw{idx};
        sta   = stationList.name{idx};
        chan  = stationList.chan{idx};
        
        
    end
end


% If station was not found, return empties
if (~flg_found)
    stLat = [];
    stLon = [];
    nw    = [];
    sta   = [];
    chan  = [];
end