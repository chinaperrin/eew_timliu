function [eqType,eqDate,t0,eqLat,eqLon,eqZ,magn,magnType,pxList] = read_phase_file(phase_file_name,o_header_only)
% Reads stp-pick file and returns eq-origin time in seconds since midnight
% as well as lists for all entries: network, station, channel, pick-time,
% epicentral distance and phase.
%
% mam 121122

% pxList: [network, channel, station, phase, pick-time]

% Units:    eqZ     [km]
%           stAlt   [m]

%phase_file_name = '../../data/wforms/brawley_2012/15198585/15198585.phase';

if (nargin < 2)
    o_header_only = 0;
end

fid    = fopen(phase_file_name);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

% End of line
eolLim = '\r?\n';

% Line
LinePattern = ['[^\r\n]*', eolLim];
thelines    = regexp(StrRay,LinePattern,'match');
nlines         = numel(thelines);

% Get header line and then discard cell entry
headerLine = thelines{1};
headerCell = textscan(headerLine,'%s');
headerCell = headerCell{1};

eqType     = headerCell{2};

eqTime     = textscan(headerCell{3},'%s','delimiter',',');
eqDate     = eqTime{1}{1};
t0         = eqTime{1}{2};

eqLat      = str2num(headerCell{4});
eqLon      = str2num(headerCell{5});
eqZ        = str2num(headerCell{6});
magn       = str2num(headerCell{7});
magnType   = headerCell{8};


if (~o_header_only)

    % Go through all lines beneath header ...
    thelines  = thelines(2:end);
    nlines    = nlines - 1;
    
    pxList.nw       = cell(nlines,1);
    pxList.sta      = cell(nlines,1);
    pxList.chan     = cell(nlines,1);
    pxList.phase    = cell(nlines,1);
    pxList.pxTime   = zeros(nlines,1);
    pxList.epiDist  = zeros(nlines,1);
    
    pxList.stLat    = zeros(nlines,1);
    pxList.stLon    = zeros(nlines,1);
    pxList.stAlt    = zeros(nlines,1);
    
    for iline = 1:nlines
        
        currentLine      = thelines{iline};
        
        pxList.nw{iline}       = currentLine(2:3);                  % NW
        pxList.sta{iline}      = strtrim(currentLine(5:10));        % Station
        pxList.chan{iline}     = currentLine(12:14);                % Band Code
        pxList.phase{iline}    = currentLine(48);                   % Phase of pick 
        pxList.pxTime(iline)   = str2double(currentLine(69:74));    % Pick time
        pxList.epiDist(iline)  = str2double(currentLine(61:66));    % Epicentral Distance
        
        pxList.stLat(iline)    = str2double(currentLine(21:27));
        pxList.stLon(iline)    = str2double(currentLine(30:38));
        pxList.stAlt(iline)    = str2double(currentLine(40:46));

    end
end
