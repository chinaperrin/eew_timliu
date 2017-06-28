function [cmt] = import_cmtList(ListFileName)

% Note: the Mw-magnitudes in fNetList saturate, e.g. Tohoku is a 8.7
%       use the JMA-magnitude instead.
%       later comment: Tohoku probably saturated because only Japanese,
%       i.e. near-source stations were used. This is not gonna be better
%       for JMA mags, right?


fid    = fopen(ListFileName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

% End of line
eolLim = '\r?\n';

% Line
LinePattern = ['[^\r\n]*', eolLim];
thelines    = regexp(StrRay,LinePattern,'match');
neq         = numel(thelines);

cmt.date  = cell(neq,1);
cmt.time  = cell(neq,1);
cmt.eqLat = zeros(neq,1);
cmt.eqLon = zeros(neq,1);
cmt.eqZ   = zeros(neq,1);
cmt.mb    = zeros(neq,1);
cmt.MS    = zeros(neq,1);

for ieq=1:neq
    
    ctLine  = thelines{ieq};
    fields  = regexp(ctLine,'\s+', 'split');
    
    cmt.date{ieq}  = fields{2};
    cmt.time{ieq}  = fields{3};
    cmt.eqLat(ieq) = str2double(fields{4});
    cmt.eqLon(ieq) = str2double(fields{5});
    cmt.eqZ(ieq)   = str2double(fields{6});
    cmt.mb(ieq)    = str2double(fields{7});
    cmt.MS(ieq)    = str2double(fields{8});
end