function out = read_noaa_buoy_realtime_files(fileFullName)

fid    = fopen(fileFullName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

eolLim      = '\r?\n';              % End of line
LinePattern = ['[^\r\n]*', eolLim]; % Line
thelines    = regexp(StrRay,LinePattern,'match');
nlines      = numel(thelines);

out.datnum = zeros(nlines-2,1);
out.wvht   = zeros(nlines-2,1);
out.swh    = zeros(nlines-2,1);
out.swp    = zeros(nlines-2,1);
out.wwh    = zeros(nlines-2,1);
out.wwp    = zeros(nlines-2,1);

for iline = 3:nlines
    
    thisLine = thelines{iline};
    fields           = regexp(thisLine,' ','split')';

    yr = str2double(fields{1});
    mt = str2double(fields{2});
    dy = str2double(fields{3});
    hr = str2double(fields{4});
    mn = str2double(fields{5});
    sc = 0;
    out.datnum(iline-2) = datenum(yr,mt,dy,hr,mn,sc);
    
    out.wvht(iline-2) = str2double(fields{7});
    out.swh (iline-2) = str2double(fields{9});
    out.swp (iline-2) = str2double(fields{10});
    out.wwh (iline-2) = str2double(fields{12});
    out.wwp (iline-2) = str2double(fields{14});
end

% Make oldest entry come first
out.datnum = flipud(out.datnum);
out.wvht   = flipud(out.wvht);
out.swh    = flipud(out.swh);
out.swp    = flipud(out.swp);
out.wwh    = flipud(out.wwh);
out.wwp    = flipud(out.wwp);