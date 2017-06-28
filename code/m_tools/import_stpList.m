function [stpList] = import_stpList(ListFileName)


fid    = fopen(ListFileName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

% End of line
eolLim = '\r?\n';

% Line
LinePattern = ['[^\r\n]*', eolLim];
thelines    = regexp(StrRay,LinePattern,'match');
neq         = numel(thelines);

stpList.id     = zeros(neq,1);
stpList.eqtype = cell(neq,1);
stpList.date   = cell(neq,1);
stpList.eqLat  = zeros(neq,1);
stpList.eqLon  = zeros(neq,1);
stpList.eqZ    = zeros(neq,1);
stpList.m      = zeros(neq,1);
stpList.mtype  = cell(neq,1);

for ieq=1:neq

      ctLine  = thelines{ieq};
      fields  = regexp(ctLine,'\s+', 'split');
      
      stpList.id(ieq)     = str2double(fields{2});
      stpList.eqtype{ieq} = fields{3};
      stpList.date{ieq}   = fields{4};
      
      stpList.eqLat(ieq)  = str2double(fields{5});
      stpList.eqLon(ieq)  = str2double(fields{6});
      stpList.eqZ(ieq)    = str2double(fields{7});

      stpList.m(ieq)      = str2double(fields{8});
      stpList.mtype{ieq}  = fields{9};
      
end