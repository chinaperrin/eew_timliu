function [station] = read_stFile_fB(stFileName)

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
station.ch     = cell(nst,1);
station.stLat  = zeros(nst,1);
station.stLon  = zeros(nst,1);

for i=1:nst
    ctLine             = thelines{i};
    station.nw{i}      = ctLine(1:2);
    station.name{i}    = ctLine(5:9);
    station.ch{i}      = ctLine(11:13);
    station.stLat(i,1) = str2num(ctLine(52:59));
    station.stLon(i,1) = str2num(ctLine(61:70));
end
station.name = strtrim(station.name);


% Clean out doubles from list
[~,ind_uniqueSt,~] =unique(station.name,'first');

station.nw    = station.nw(ind_uniqueSt);
station.name  = station.name(ind_uniqueSt);
station.ch    = station.ch(ind_uniqueSt);
station.stLat = station.stLat(ind_uniqueSt);
station.stLon = station.stLon(ind_uniqueSt);

station.fullname = strcat(station.nw,'.',station.name);

fprintf(1,[' ',num2str(numel(ind_uniqueSt)), ' unique stations found\n'])

