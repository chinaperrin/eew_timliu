function [sraw,meta] = read_ascii_wform_wenchuan(traceFullName)

nhdr = 16;
% Reads Wenchuan Ascii Files that I have gotten from Chaoyong Peng via
% researchGate, and files from tar-file from Lingling
%
% menandrin@gmail.com, 161122

fid    = fopen(traceFullName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

eolLim      = '\r?\n';                              % End of line
LinePattern = ['[^\r\n]*', eolLim];                 % Line
thelines    = regexp(StrRay,LinePattern,'match');

meta.station.unit = 'cm/s/s';


tmp=regexp(thelines{6},' ','split');
meta.station.name = strtrim(tmp{end});

tmp=regexp(thelines{7},' ','split');
meta.station.ground   = strtrim(tmp{end});

tmp=regexp(thelines{8},' ','split');
meta.station.itype   = strtrim(tmp{end});

tmp=regexp(thelines{10},' ','split');
meta.station.component   = strtrim(tmp{end});

tmp=regexp(thelines{12},' ','split');
meta.record.npoints = str2double(tmp{5});
meta.station.sr     = 1/str2double(tmp{12});

tmp=regexp(thelines{13},' ','split');
meta.record.maxamp   = str2double(tmp{4});
meta.record.tmax     = str2double(tmp{10});
meta.record.duration = str2double(tmp{17});

tmp=regexp(thelines{14},' ','split');
meta.record.tpre   = str2double(tmp{4});



% Read waveform amplitudes
fullLines = thelines(nhdr+1:end-1);
nlines    = numel(fullLines);
raw       = str2num(cell2mat(fullLines));
ncol      = size(raw,2);
sraw      = reshape(raw',nlines*ncol,1);
 

% Add last line
lastline = thelines(end);
lastline = str2num(cell2mat(lastline));
sraw     = [sraw;lastline'];