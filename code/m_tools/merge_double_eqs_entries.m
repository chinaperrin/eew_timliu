function eqs = merge_double_eqs_entries(eqs,eqName1,eqName2)

% Add traceId's etc from eqName2 to eqName1, and then delete eqName2 from
% list
i1 = find(strcmp(eqs.name,eqName1));
i2 = find(strcmp(eqs.name,eqName2));

if ~isempty(i1) &~isempty(i2)
    
    % Add entries from second entry to those of first one
    neqs            = numel(eqs.m);
    eqs.traceId{i1} = [eqs.traceId{i1}; eqs.traceId{i2}];
    eqs.tpx    {i1} = [eqs.tpx{i1}    ; eqs.tpx{i2}];
    eqs.lt     {i1} = [eqs.lt{i1}     ; eqs.lt{i2}];
    
    % Remove second entry from list
    aint_i2     = setdiff(1:neqs,i2);
    eqs.eventId = eqs.eventId(aint_i2);
    eqs.traceId = eqs.traceId(aint_i2);
    eqs.tpx     = eqs.tpx    (aint_i2);
    eqs.lt      = eqs.lt     (aint_i2);
    eqs.date    = eqs.date   (aint_i2);
    eqs.t0      = eqs.t0     (aint_i2);
    eqs.name    = eqs.name   (aint_i2);
    eqs.lat     = eqs.lat    (aint_i2);
    eqs.lon     = eqs.lon    (aint_i2);
    eqs.z       = eqs.z      (aint_i2);
    eqs.m       = eqs.m      (aint_i2);
    
else
    fprintf(1,'Only one or neither entry is contained in eqs \n\t--> Proceeding without changing eqs\n')
end