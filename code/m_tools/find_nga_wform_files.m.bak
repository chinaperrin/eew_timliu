function [rawFullName,targetName] = find_nga_wform_files(rawName,wformDir,o_verbose)

if (nargin<3); o_verbose=0; end

nc = 3;   % Default value, in case no orntCode can be identified
np = 0;   % Nr. of times a horizontal name has been processed (to detect multiple processings = error)

% 1. Preprocess name: trim spaces, replace backslashes and remove extensions
tmp      = regexp(rawName,'\/','split');
dirname  = tmp{1};
fname    = strrep(tmp{2}, '.at2','');
fname    = strrep(fname,  '.AT2','');

% 2. Detect orientation code and baseName (= name without orientation code)
% Vertical    orientation codes:      -UP. -V1. DN. DWN. --Z. --V. -V. 
% NS/SN       orientation codes:      -N. NS. -NS.
% EW/WE           "         "  :      -WE. -E. 
% other horizontal    "   "    :      [0-9][0-9][0-9]. -TR -LN -T1 -L1 --X --Y


% a. vertical    -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
z_orntCodes = {'-UP';'-V1';'DN';'DWN';'--Z';'--V';'-V'};
z_match     = regexp(fname(end-2:end),z_orntCodes);
z_idx       = find(cellfun(@(x) ~isempty(x), z_match));
if (~isempty(z_idx))
    code     = z_orntCodes{z_idx};
    orntCode = 'V';
    nc       = numel(code);
    np       = np+1;
end


% b. horizontal     -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
hn_orntCodes   = {'-NS';'-SN';'NS';'SN';'-N';'-S'};         % e/w orientations
he_orntCodes   = {'-EW';'-WE';'EW';'WE';'-W';'-E'};         % n/s orientations
ho_orntCodes   = {'-TR';'-LN';'-T1';'-L1';'--X';'--Y'};     % other (ambiguous) horizontal orientations
ha_orntCodes   = {'\d\d\d'};                                  % Code in the form of an azimuth (e.g. '180' or '090')

hn_match     = regexp(fname(end-2:end),hn_orntCodes);
he_match     = regexp(fname(end-2:end),he_orntCodes);
ho_match     = regexp(fname(end-2:end),ho_orntCodes);
ha_match     = regexp(fname(end-2:end),ha_orntCodes);

hn_idx       = find(cellfun(@(x) ~isempty(x), hn_match));
he_idx       = find(cellfun(@(x) ~isempty(x), he_match));
ho_idx       = find(cellfun(@(x) ~isempty(x), ho_match));
ha_idx       = find(cellfun(@(x) ~isempty(x), ha_match));

if (~isempty(hn_idx))
    code     = hn_orntCodes{hn_idx};
    orntCode = '999';
    nc       = numel(code);
    np       = np+1;
end
if (~isempty(he_idx))
    code     = he_orntCodes{he_idx};
    orntCode = '888';
    nc       = numel(code);
    np       = np+1;
end
if (~isempty(ho_idx))
    code     = ho_orntCodes{ho_idx};
    strng    = regexpi(code,'[a-z]+','match');
    orntCode = strcat(strng{1},'_777');
    nc       = numel(code);
    np       = np+1;
end
if (~isempty(ha_idx))
    orntCode = fname(end-2:end);
    nc       = 3;
    np       = np+1;
end

if (np>1) fprintf(1,'8UNG: fname has matched more than one pattern. WHATNOW?\n'); pause; end

% 3. Find file by searching for baseName
baseName = fname(1:end-nc);

if (o_verbose)
    fprintf(1,['   Searching for ',rawName,' by the name of ',fname,'* in dir ',dirname,'\n\t']);
end

[~,result]  = unix(['find ',wformDir, ' -name ', fname,'*']);
eolLim      = '\r?\n';              % End of line
LinePattern = ['[^\r\n]*', eolLim]; % Line
resultlines = regexp(result,LinePattern,'match');

% Find the one with the right directory name
match    = regexp(resultlines,dirname);
idxMatch = find(cellfun(@(x) ~isempty(x), match));

if (numel(idxMatch)==1)
    rawFullName = strtrim(resultlines{idxMatch});
    if (o_verbose); fprintf(1,['One match: ',rawFullName,'\n']); end
   
    % Create unique filename
    targetName = strcat([dirname,'_',baseName,'.',orntCode,'.AT2']);

elseif (numel(idxMatch)>1)
    fprintf(1,['Multiple matches: ',resultlines{idxMatch},' \n\n\t WHATNOW?'])
    pause
    
elseif (numel(idxMatch)==0)
    if (o_verbose); fprintf(1,['No matches: ',result,'\n']); end
    rawFullName = [];
    targetName  = [];
end