function list = read_list_with_gaps_kikuchi(listFullName)

fid    = fopen(listFullName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

eolLim      = '\r?\n';
LinePattern = ['[^\r\n]*', eolLim];
thelines    = regexp(StrRay,LinePattern,'match');
thelines    = thelines(3:end);

ndat      = numel(thelines);
list.date = zeros(ndat,1);
list.Mw   = zeros(ndat,1);
list.M0   = zeros(ndat,1);
list.T    = zeros(ndat,1);
list.z    = zeros(ndat,1);

for idat=1:ndat
    
    ctLine = thelines{idat};
    
    list.date(idat) = str2double(ctLine(1:8));
    list.Mw  (idat) = str2double(ctLine(46:48));
    list.M0  (idat) = str2double(ctLine(50:56));
    list.T   (idat) = str2double(ctLine(60:64));
    list.z   (idat) = str2double(ctLine(38:40));
end