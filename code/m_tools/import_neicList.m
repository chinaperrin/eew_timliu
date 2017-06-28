function [neicList] = import_neicList(ListFileName)


fid    = fopen(ListFileName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

% End of line
eolLim = '\r?\n';

% Line
LinePattern = ['[^\r\n]*', eolLim];
thelines    = regexp(StrRay,LinePattern,'match');
neq         = numel(thelines);

neicList.date   = cell(neq,1);
neicList.time   = cell(neq,1);
neicList.eqLat  = zeros(neq,1);
neicList.eqLon  = zeros(neq,1);
neicList.eqZ    = zeros(neq,1);
neicList.m      = zeros(neq,1);
neicList.mtype  = cell(neq,1);
neicList.sc     = cell(neq,1);

for ieq=1:neq

      ctLine  = thelines{ieq};
      fields  = regexp(ctLine,'\s+', 'split');
      
      neicList.date{ieq}   = fields{1};
      neicList.time{ieq}   = fields{2};
      
      neicList.eqLat(ieq)  = str2double(fields{3});
      neicList.eqLon(ieq)  = str2double(fields{4});
      neicList.eqZ(ieq)    = str2double(fields{5});

      neicList.m(ieq)      = str2double(fields{6});
      neicList.mtype{ieq}  = fields{7};
      neicList.sc{ieq}     = fields{8};
end