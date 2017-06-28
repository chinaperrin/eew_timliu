function [evntList] = import_fnetList(ListFileName)

% Note: the Mw-magnitudes in fNetList saturate, e.g. Tohoku is a 8.7
%       use the JMA-magnitude instead.


fid    = fopen(ListFileName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

% End of line
eolLim = '\r?\n';

% Line
LinePattern = ['[^\r\n]*', eolLim];
thelines    = regexp(StrRay,LinePattern,'match');
neq         = numel(thelines);

evntList.date   = cell(neq,1);
evntList.eqLat  = zeros(neq,1);
evntList.eqLon  = zeros(neq,1);
evntList.jmaZ   = zeros(neq,1);
evntList.jmaM   = zeros(neq,1);
evntList.strike = cell(neq,1);
evntList.dip    = cell(neq,1);
evntList.rake   = cell(neq,1);
evntList.M0     = zeros(neq,1);
evntList.Mw     = zeros(neq,1);
evntList.MTZ    = zeros(neq,1,'single');

for ieq=1:neq

      ctLine  = thelines{ieq};
      fields  = regexp(ctLine,'\s+', 'split');
      
      evntList.date{ieq}   = fields{1};
      evntList.eqLat(ieq)  = str2double(fields{2});
      evntList.eqLon(ieq)  = str2double(fields{3});
      evntList.jmaZ(ieq)   = str2double(fields{4});
      evntList.jmaM(ieq)   = str2double(fields{5});

      evntList.strike{ieq} = fields{7};
      evntList.dip{ieq}    = fields{8};
      evntList.rake{ieq}   = fields{9};
      
      evntList.M0(ieq)     = str2double(fields{10});
      evntList.Mw(ieq)     = str2double(fields{12});
      evntList.MTZ(ieq)    = str2double(fields{11});
end