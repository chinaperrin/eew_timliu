function [epiDist] = get_epiDist(nw,sta,eqLat,eqLon)

global stationList missingStationList

% Find station in list
stFullname = strcat(nw,'.',sta);
stIdx      = find(strcmp(stFullname,stationList.fullname));
%stIdx      = find(strcmp(stFullname,stationList.fullname),1,'first');

if (numel(stIdx)>1)
    fprintf(1,['8UNG:  Station ',nw,'.',sta,': More than one match in stationList. 1st one used.\n'])
    stIdx = stIdx(1);
    multipleMatchList = [multipleMatchList; strcat(nw,'.',sta)];
end


if (isempty(stIdx))
    
    % fprintf(1,['8UNG:  Station ',nw,'.',sta,' not in list, no epicentral distance computed'])
    epiDist            = 999999;
    missingStationList = [missingStationList; strcat(nw,'.',sta)];
        
else
    
	% Compute epicenter [m]
    stLat   = stationList.stLat(stIdx);
    stLon   = stationList.stLon(stIdx);
    epiDist = vdist(eqLat,eqLon,stLat,stLon);
end