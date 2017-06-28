function [TraceList] = delete_blackListed_records(TraceList,blackListFullName,opts)


load(blackListFullName)

nbk     = numel(blackList.m);
iDelete = [];

for ibk = 1:nbk
    
    % Find traceIdx and indices of other comps of same record
    fullName         = blackList.fullName{ibk};
    [idxZ,idxE,idxN] = find_corec_idx(fullName,TraceList,0);
    
    ntrfound = sum(cellfun(@(x) ~isempty(x), {idxZ; idxE; idxN}));
    if ntrfound==3
    
        % Save traceList indices
        iDelete = [iDelete; idxZ; idxE; idxN];
        m = TraceList.m(idxZ);
        r = TraceList.hypDist(idxZ);

        fprintf(1,sprintf('\n%i/%i: Removing 3 traces from TraceList (m%3.1f@%5.1fkm); reason: %s',ibk,nbk,m,r,blackList.comment{ibk}))
        TraceList.fullName([idxZ; idxE; idxN])
        
        
        % Optionally plot traces
        if opts.plot_blackListedTraces
            
            TraceList.printSingleTraceSummary(idxZ)
            
            % Read vertical trace and plot raw waveform
            [sz,meta] = read_any_trace(TraceList.fullName{idxZ},TraceList,opts);
            [se,~   ] = read_any_trace(TraceList.fullName{idxZ},TraceList,opts);
            [sn,~   ] = read_any_trace(TraceList.fullName{idxZ},TraceList,opts);
            t         = meta.t;
            ppxIdx    = meta.ppxIdx;
            tppx      = t(ppxIdx);
            
            % Plot with pick from traceList, vel-pick and acc-pick
            sig.szraw  = sz.raw; % Raw, only gain corrected
            sig.seraw  = se.raw; % Raw, only gain corrected
            sig.snraw  = sn.raw; % Raw, only gain corrected
            allPx.tppx = tppx;
            plot_wform_components(sig,t,allPx,[tppx-.2 tppx+.05],[])
        end
        
    else
        fprintf(1,'At least one component of black-listed entry not found in traceList. \nNot deleting these records from TraceList.\n')
    end
    1+1;
end

% Remove all identified records from TraceList
idxKeep   = setdiff(1:numel(TraceList.m)',iDelete)';
TraceList = TraceList.selectSubList(idxKeep);
