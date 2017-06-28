function [doubleIdx] = find_double_records(TraceList)

% Find indices of strong motion records that have a colocated broadband
% record. Assumes clipped BB waveforms have already been removed from
% TraceList, i.e. if both are still in traceList the SM record can be
% romved

o_debug = false;

%[TraceList,SkipList,fc,snpLength] = import_fbOut({'/scratch/memeier/data/socal/scsn_900101_011231/out/i32/acausal/'});

ntr       = numel(TraceList.eq.m);
itr       = 1;
doubleIdx = [];

while itr<ntr
    
    print_iteration_numbers(itr,ntr,'hundreds')
    
    % Find all records from same site
    fullName   = TraceList.fullName{itr};
    [iZ,iE,iN] = find_corec_idx(fullName,TraceList,1);
    iAll       = sort([iZ;iE;iN]);
    
    % Identify strong motion stations
    isSM  = regexp(TraceList.fullName(iAll),'H[LN][ENZ]','match');
    idxSM = find(cellfun(@(x) ~isempty(x), isSM));
    nsm   = numel(idxSM);

    % Identify broad band stations
    isBB  = regexp(TraceList.fullName(iAll),'HH[ENZ]','match');
    idxBB = find(cellfun(@(x) ~isempty(x), isBB));
    nbb   = numel(idxBB);

    % If there are 3 of either instrument type, skip the SM ones. (The fact 
    % that the BB record is still in the traceList at this point implies
    % that itr is not clipped
    if ( (nsm==3) && (nbb==3) )
        doubleIdx = [doubleIdx; iAll(idxSM)];
        
        % Read and plot waveforms
        if o_debug
            
            % Read z-components
            i1 = iZ(1);
            [zS1,meta] = read_any_trace(TraceList.fullName{i1},TraceList,1);
            t1         = meta.t(1:numel(zS1.vel));
            plot_wform(zS1.vel,t1,TraceList.tppx(i1),TraceList.tspx(i1),194);

            i2 = iZ(2);
            [zS2,meta] = read_any_trace(TraceList.fullName{i2},TraceList,1);
            t2         = meta.t(1:numel(zS2.vel));
            hold on; plot(t2,zS2.vel,':r')
            
        end
    end
    
    itr = max(iAll)+1;
end