function [station] = read_FDSN_stationFile(fileFullName)

fid    = fopen(fileFullName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

eolLim      = '\r?\n';   % End of line
LinePattern = ['[^\r\n]*', eolLim];
thelines    = regexp(StrRay,LinePattern,'match');

if ~isempty(thelines);
    fields             = regexp(thelines{2},'\|','split');
    station.network    = fields{1};
    station.name       = fields{2};
    station.location   = fields{3};
    station.channel    = fields{4};
    station.lat        = str2double(fields{5});
    station.lon        = str2double(fields{6});
    station.elevation  = str2double(fields{7});
    station.depth      = str2double(fields{8});
    station.azimuth    = str2double(fields{9});
    station.dip        = str2double(fields{10});
    station.instrument = fields{11};
    station.scale      = str2double(fields{12});
    station.scaleFreq  = str2double(fields{13});
    station.scaleUnit  = fields{14};
    
else
    station = [];
end