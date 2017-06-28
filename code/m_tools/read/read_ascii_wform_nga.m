function [sraw,meta] = read_ascii_wform_nga(traceFullName,o_hdr)

% Reads NGA West1 ascii waveform files
% menandrin@gmail.com, 130912

%tic
% if (nargin<3)
%     o_hdr = 0;
% end

% nhdr  = 4;
% o_hdr = 1;
% traceFullName = '~/Downloads/sucker/peer.berkeley.edu/nga_files/ath/LANDERS/GR2-UP.AT2';

nhdr = 4;

fid    = fopen(traceFullName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

eolLim      = '\r?\n';                  % End of line
LinePattern = ['[^\r\n]*', eolLim];     % Line pattern
thelines    = regexp(StrRay,LinePattern,'match');


% Read header lines
if (~o_hdr)
    meta = [];
else
    
    hdrLines  = thelines(1:nhdr);
    
    tmp       = regexp(hdrLines{4},'[\ ]*','split');
    dt        = str2num(tmp{2});
    meta.sr   = 1/dt;
    
    tmp       = regexp(hdrLines{3},'\ ','split');
    meta.unit = tmp{end}(1);
end



% Read waveform amplitudes
lastline    = thelines(end);
lastline    = str2num(cell2mat(lastline));
thelines    = thelines(nhdr+1:end-1);
nlines      = numel(thelines);

raw  = str2num(cell2mat(thelines));
ncol = size(raw,2);
sraw = reshape(raw',nlines*ncol,1);
sraw = [sraw;lastline'];