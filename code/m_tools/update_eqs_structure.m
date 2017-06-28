function eqs = compile_eqs_structure(trList,skipList)
% Takes existing <eqs> structure that had been created with the function 
% compile_eqs_structure.m and updates
global vp

thousands = linspace(1e3,1e6,1e3);

% A. Assign eqIdx to all entries of trList ................................
ntr          = numel(trList.m);
trList.eqIdx = zeros(ntr,1,'int32');
id           = 0;

for itr = 1:ntr
    if (trList.eqIdx(itr)==0)
        fullName        = trList.fullName{itr};
        [ids_thisEvent] = find_all_event_traces(fullName,trList);
        if ~isempty(skipList)
            [sids_thisEvent] = find_all_event_traces(fullName,skipList);
        end
            
        % Sometimes find_all_event_traces.m misclassifies traces. Check
        % if all associated traces have same magnitudes and if not,
        % assign them to different events
        if ~isempty(skipList) mVals = unique([trList.m(ids_thisEvent); skipList.m(sids_thisEvent)]);
        else                  mVals = unique(trList.m(ids_thisEvent));
        end
        
        for im = 1:numel(mVals)
            
            id = id+1;
            
            % TraceList
            itmp              = find(mVals(im) == trList.m(ids_thisEvent));
            idx               = ids_thisEvent(itmp);
            trList.eqIdx(idx) = id;
            
            % SkipList
            if ~isempty(skipList); itmps                = find(mVals(im) == skipList.m(sids_thisEvent));
                                   idxs                 = sids_thisEvent(itmps);
                                   skipList.eqIdx(idxs) = id;
                 fprintf(1,sprintf('id = %i --- %i  records (+%i skipped) --- trList-index = %i/%i\n',id,numel(idx),numel(idxs),itr,ntr))
            else fprintf(1,sprintf('id = %i --- %i  records --- trList-index = %i/%i\n',id,numel(idx),itr,ntr))
            end
        end
    end
end


% B. Unique list of eqIds .................................................
% <eqs.eventId> = list of all event-idndices considered for inference
eqs.eventId = unique(trList.eqIdx);
neq         = numel(eqs.eventId);


% C. Save trace indices for all traces in each eqs.eventId ................
%    sorted wrt/ distance, along with the corresp. theretical p-arrival times
eqs.traceId = cell(neq,1);
eqs.tpx     = cell(neq,1);
eqs.lt      = cell(neq,1);
eqs.date    = cell(neq,1);
eqs.t0      = cell(neq,1);
%eqs.mType   = cell(neq,1);
eqs.lat     = zeros(neq,1);
eqs.lon     = zeros(neq,1);
eqs.z       = zeros(neq,1);

for ieq = 1:neq
    if ismember(ieq,thousands); fprintf(1,[num2str(ieq),'/',num2str(neq),'\n']);end
    
    trIdx            = find(eqs.eventId(ieq)==trList.eqIdx);      % traceList index of all traces from target event
    R                = trList.hypDist(trIdx);
    [Rstd,srtIdx]    = sort(R);
    tp               = Rstd/vp;
    
    eqs.traceId{ieq} = trIdx(srtIdx);
    eqs.tpx{ieq}     = tp;
    eqs.date{ieq}    = trList.eqDate{trIdx(1)};
    eqs.t0{ieq}      = trList.t0    {trIdx(1)};
    %eqs.mType{ieq}   = trList.mType {trIdx(1)};
    eqs.lat(ieq)     = trList.eqLat (trIdx(1));
    eqs.lon(ieq)     = trList.eqLon (trIdx(1));
    eqs.z(ieq)       = trList.eqZ   (trIdx(1));
    
    if ~isempty(skipList)
        skipIdx               = find(eqs.eventId(ieq)==skipList.eqIdx);
        Rskip                 = skipList.hypDist(skipIdx);
        [RskipStd,skipSrtIdx] = sort(Rskip);
        tpskip                = RskipStd/vp;
        
        eqs.skipId {ieq} = skipIdx(skipSrtIdx);
        eqs.Rskip{ieq}   = RskipStd;
        eqs.tpxSkip{ieq} = tpskip;
    end

end
eqs.eqIdx = trList.eqIdx;
eqs.m     = trList.m(cellfun(@(x) x(1), eqs.traceId));