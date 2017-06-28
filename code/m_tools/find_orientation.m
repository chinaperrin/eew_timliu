function [orientation,instrument] = find_orientation(GL)


ntr            = numel(GL.fullName);
traceNames_raw = GL.fullName;
dsNames_raw    = GL.dataSetName;

orientation = cell(ntr,1);
instrument  = cell(ntr,1);

for itr = 1:ntr

    trName        = traceNames_raw{itr};

    % JAPAN
    if ( (strcmp(dsNames_raw(itr),'k_kik_Fnet')) || (strcmp(dsNames_raw(itr),'k_kik_noFnet')) )

        % So far using only strong motion data from Japan
        instrument{itr}  = 'L';

        % Sensor orientation
        ptIdx = regexp(trName,'\.');
        o_tmp = trName(ptIdx+1:ptIdx+2);

        if     (strcmp(o_tmp,'UD')); orientation{itr} = 'Z';
        elseif (strcmp(o_tmp,'NS')); orientation{itr} = 'N';
        elseif (strcmp(o_tmp,'EW')); orientation{itr} = 'E';
        else   fprintf(1,'orientation NOT FOUND ... what now?\n'); pause;
        end


    % CALIFORNIA
    elseif ( (strcmp(dsNames_raw(itr),'scsnPx')) || (strcmp(dsNames_raw(itr),'scsnTs')) )
        sacIdx = regexp(trName,'.sac');

        instrument{itr}  = trName(sacIdx-2);
        orientation{itr} = trName(sacIdx-1);


    % NOT FOUND     
    else
        fprintf(1,'NOT FOUND ... what now?\n')
        pause
    end
end