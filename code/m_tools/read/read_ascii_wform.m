function [sraw] = read_ascii_wform(traceFullName,nhdr,ncol)

fid    = fopen(traceFullName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

% End of line
eolLim = '\r?\n';

% Line
LinePattern = ['[^\r\n]*', eolLim];
thelines    = regexp(StrRay,LinePattern,'match');
thelines    = thelines(nhdr+1:end);
nlines      = numel(thelines);

raw  = str2num(cell2mat(thelines));
sraw = reshape(raw',nlines*ncol,1);