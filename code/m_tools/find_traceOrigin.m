function [traceOrigin] = find_traceOrigin(traceFullName)

% Find out what sort of waveform it is
ptIdx     = regexp(traceFullName,'\.');
extension = traceFullName(ptIdx(end):end);

jp_extensions = {'.EW','.EW1','.EW2','.NS','.NS1','.NS2','.UD','.UD1','.UD2'};

iwenchuan1 = regexp(traceFullName,'_ORIG_','ONCE');
iwenchuan2 = regexp(traceFullName,'/china/','ONCE');

% JAPAN - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ismember(extension,jp_extensions)
    traceOrigin = 'japan';
    
elseif strcmp(extension,'.sac')
    traceOrigin = 'cali';

elseif strcmp(extension,'.mseed')
    traceOrigin = 'nocal';

elseif strcmp(extension,'.AT2')
    traceOrigin = 'nga';

elseif ~isempty(iwenchuan1) &&~isempty(iwenchuan2) 
    traceOrigin = 'wenchuan';
    
else
    traceOrigin = [];
    fprintf('Not sure where this trace is coming from ... What now?\n')
end