function [station] = import_stationlist_fB(stFileName)

fid    = fopen(stFileName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

% End of line
eolLim = '\r?\n';

% Line
LinePattern = ['[^\r\n]*', eolLim];
thelines    = regexp(StrRay,LinePattern,'match');
nst         = numel(thelines);

station.nw      = cell(nst,1);
station.name    = cell(nst,1);
station.lat   = zeros(nst,1);
station.lon   = zeros(nst,1);
station.ttt     = cell(nst,1);

for i=1:nst

    ctLine             = thelines{i};
    
    station.name{i}    = ctLine(1:5);
    station.nw{i}      = ctLine(7:8);
    
    % Station oordinates
    tmp_lat1           = ctLine(14:15);
    tmp_lat2           = ctLine(17:23);
    station.lat(i)   = str2num(tmp_lat1) + str2num(tmp_lat2)/60;
    
    tmp_lon1           = ctLine(25:27);
    tmp_lon2           = ctLine(29:35);
    station.lon(i)   = (str2num(tmp_lon1) + str2num(tmp_lon2)/60)*-1;
    
end

station.name = strtrim(station.name);

% Clean out doubles from list (10 stations from BK-network are also in cinnabar-list)
%[~,ind_uniqueSt,~] =unique(station.name,'first');
[~,sorted_idx,~] = unique(station.name,'first');
ind_uniqueSt     = sort(sorted_idx);

station.nw    = station.nw(ind_uniqueSt);
station.name  = station.name(ind_uniqueSt);
station.lat   = station.lat(ind_uniqueSt);
station.lon   = station.lon(ind_uniqueSt);
station.ttt   = station.ttt(ind_uniqueSt);

nst           = numel(station.lat(:,1));
fprintf(1,[' ... ',num2str(numel(ind_uniqueSt)), ' out of ',num2str(nst),' stations in California station list are unique\n'])
