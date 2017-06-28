function [match] = find_colocated_SM_BB_records(trList)

ntr = numel(trList.m);
hasBeenProcessed = false(ntr,1);

thousands = linspace(1e3,1e6,1e6/1e3);

dmax_stLat = 0;
dmax_stLon = 0; 
dmax_date  = 0;
dmax_m     = 0;
dmax_r     = 1;

eqLat  = trList.eqLat;
eqLon  = trList.eqLon;
stLat  = trList.stationLat;
stLon  = trList.stationLon;
m      = trList.m;
r      = trList.hypDist;
eqDate = datenum(trList.eqDate);

match.bbIdx = cell(0);
match.smIdx = cell(0);

for itr = 1:ntr
    
    if ismember(itr,thousands); fprintf(1,sprintf('%i / %i ..\n',itr,ntr)); end
    
    % If trace hasn't already been assigned
    if ~hasBeenProcessed(itr)
        dstLat  = stLat(itr)  -  stLat;
        dstLon  = stLon(itr)  -  stLon;
        deqDate = eqDate(itr) - eqDate;
        dm      = m(itr)      - m;
        dr      = r(itr)      - r;
        
        % Indices of traces with same date and station coordinates
        idx = find(abs(dstLat)<=dmax_stLat &abs(dstLon)<=dmax_stLon &abs(deqDate)<=dmax_date &abs(dm)<=dmax_m &abs(dr)<=dmax_r);
        idx = setdiff(idx,itr);
        
        % If there are potential matches ...
        if ~isempty(idx)
            
            % If investigated trace is a BB record ...
            if strcmp(trList.instrCode(itr),'H')

                % ... add all found indices to match-list that are SM records
                idx = idx(~strcmp(trList.instrCode(idx),'H'));
                if ~isempty(idx)

                    %fprintf(1,'\n\nMatch!')
                    %trList.printSingleTraceSummary(itr);
                    %for i=1:numel(idx)
                    %    trList.printSingleTraceSummary(idx(i));
                    %end
            
                    match.bbIdx           = [match.bbIdx; itr];
                    match.smIdx           = [match.smIdx; idx];
                    hasBeenProcessed(idx) = true;
                end
                
                
            % If investigated trace is a SM record ...
            else
                
                % ... add all found indices to match-list that are BB records
                idx = idx(strcmp(trList.instrCode(idx),'H'));
                if ~isempty(idx)
                    
                    %fprintf(1,'\n\nMatch!\n')
                    %trList.printSingleTraceSummary(itr);
                    %for i=1:numel(idx)
                    %    trList.printSingleTraceSummary(idx(i));
                    %end
                    
                    match.smIdx           = [match.smIdx; itr];
                    match.bbIdx           = [match.bbIdx; idx];
                    hasBeenProcessed(idx) = true;
                end
            end
            
        end
        
    end
    hasBeenProcessed(itr) = true;
    
end