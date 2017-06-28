function quakeName = get_quake_name_oldTraceLists(traceFullName);
% Pre-i38 versions of traceLists did not have a field eqName. For those,
% use this function to cunstruct an event name, depending on the data
% source.


traceOrigin = find_traceOrigin(traceFullName);

if     strcmp(traceOrigin,'japan');    quakeName = 'Japan';
elseif strcmp(traceOrigin,'cali');     quakeName = 'Southern California';
elseif strcmp(traceOrigin,'wenchuan'); quakeName = 'China';
elseif strcmp(traceOrigin,'nga');      quakeName = 'Taiwan';
end