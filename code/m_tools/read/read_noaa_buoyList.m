function buoyList = read_noaa_buoyList(fileFullName)

nhdrLines = 7;

fid    = fopen(fileFullName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

eolLim      = '\r?\n';              % End of line
LinePattern = ['[^\r\n]*', eolLim]; % Line
thelines    = regexp(StrRay,LinePattern,'match');
nlines      = numel(thelines);

buoyList.name = cell(nlines-nhdrLines,1);
buoyList.lat  = zeros(nlines-nhdrLines,1);
buoyList.lon  = zeros(nlines-nhdrLines,1);
buoyList.z    = zeros(nlines-nhdrLines,1);

for iline = nhdrLines+1:nlines
    
    iout = iline-nhdrLines;
    
    thisLine      = thelines{iline};
    thisCleanLine = regexprep(thisLine,' +',' ');
    fields        = regexp(thisCleanLine,' ','split')';

    buoyList.name{iout} = fields{2};
    buoyList.lat (iout) = str2double(fields{4}) + str2double(fields{5})/60 + str2double(fields{6})/3600;
    buoyList.lon (iout) = -(str2double(fields{8}) + str2double(fields{9})/60 + str2double(fields{10})/3600);
    buoyList.z   (iout) = str2double(fields{12});
end