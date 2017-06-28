function [has_pPx,allPx] = get_pPick(pxList,networkName,stationName)

% If there is a broadband pick, prefer it, if there is not, take any other
% p-pick on a vertical component.

lgc_nw     = single(strcmp(networkName,pxList.nw));
lgc_sta    = single(strcmp(stationName,pxList.sta));
lgc_pPhase = single(strcmp('P',        pxList.phase));
lgc_bbPx   = single(strcmp('HHZ',      pxList.chan));   % Vertical broadband pick

tmp        = regexp(pxList.chan,'Z');                    
lgc_zPx    = cellfun(@(x) ~isempty(x), tmp);            % Any vertical pick



bbPxIdx    = find(lgc_nw.*lgc_sta.*lgc_pPhase.*lgc_bbPx);   % Vertical broadband p-pick
otherPxIdx = find(lgc_nw.*lgc_sta.*lgc_pPhase.*lgc_zPx);    % Any vertical p-pick


% If there is at least one broadband p-pick (sometimes there are more than
% one p-pick on the same trace ...)
if (~isempty(bbPxIdx))
    
    ppxIdx  = bbPxIdx(1);
    has_pPx = true;
    
% If there is no vertical p-pick on a broadband channel, choose any other vertical p-pick
elseif ( (isempty(bbPxIdx)) && (~isempty(otherPxIdx)) )
    
    ppxIdx  = otherPxIdx(1);
    has_pPx = true;
    
% If there is no vertical p-pick at all    
else
    
    ppxIdx  = [];
    has_pPx = false;
end


% Write pick information into <allPx>
if (~isempty(ppxIdx))

    allPx.nw      = pxList.nw(ppxIdx);
    allPx.sta     = pxList.sta(ppxIdx);
    allPx.chan    = pxList.chan(ppxIdx);
    allPx.phase   = pxList.phase(ppxIdx);
    allPx.tppx    = pxList.pxTime(ppxIdx);
    allPx.epiDist = pxList.epiDist(ppxIdx);
    
    allPx.stLat = pxList.stLat(ppxIdx);
    allPx.stLon = pxList.stLon(ppxIdx);
    allPx.stAlt = pxList.stAlt(ppxIdx);
    
else
    allPx   = [];
end