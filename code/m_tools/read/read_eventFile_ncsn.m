function evt = read_eventFile_ncsn(fileFullName)

fid    = fopen(fileFullName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

eolLim      = '\r?\n';              % End of line
LinePattern = ['[^\r\n]*', eolLim]; % Line
thelines    = regexp(StrRay,LinePattern,'match');
nlines      = numel(thelines);

evt.date  = cell (nlines-1,1);
evt.t0    = cell (nlines-1,1);
evt.id    = cell (nlines-1,1);
evt.lat   = zeros(nlines-1,1);
evt.lon   = zeros(nlines-1,1);
evt.depth = zeros(nlines-1,1);
evt.mtype = cell (nlines-1,1);
evt.m     = zeros(nlines-1,1);

for iline = 1:nlines-1
    
    thisLine = thelines{iline+1};
    fields           = regexp(thisLine,'\|','split');
    
    tmp              = regexp(fields{2},'T','split');
    evt.date {iline} = tmp{1};
    evt.t0   {iline} = tmp{2}; 
    
    evt.id   {iline} = fields{1};
    evt.lat  (iline) = str2double(fields{3});
    evt.lon  (iline) = str2double(fields{4});
    evt.depth(iline) = str2double(fields{5});
    evt.mtype{iline} = fields{10};
    evt.m    (iline) = str2double(fields{11});
end