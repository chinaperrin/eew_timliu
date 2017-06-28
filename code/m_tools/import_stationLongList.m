function [station] = import_stationLongList(stFileName)

fid    = fopen(stFileName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

% End of line
eolLim = '\r?\n';

% Line
LinePattern = ['[^\r\n]*', eolLim];
thelines    = regexp(StrRay,LinePattern,'match');
nst         = numel(thelines);

station.nw     = cell(nst,1);
station.name   = cell(nst,1);
station.chan   = cell(nst,1);
station.lat    = zeros(nst,1);
station.lon    = zeros(nst,1);
station.elev   = zeros(nst,1);

for i=1:nst

    ctLine             = thelines{i};
    fields             = regexp(ctLine,'\s+', 'split');
    
    station.nw{i}      = fields{1};
    station.name{i}    = fields{2};
    station.chan{i}    = fields{3};
    
    station.lat(i)   = str2double(ctLine(52:59));
    station.lon(i)   = str2double(ctLine(61:70));
    station.elev(i)  = str2double(ctLine(72:76));
    
    1+1;
end

fprintf(1,[' ... coordinates of ',num2str(nst), ' stations in California have been read\n'])
