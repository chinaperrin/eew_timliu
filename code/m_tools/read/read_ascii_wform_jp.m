function [sraw,meta] = read_ascii_wform_jp(traceFullName,nhdr,o_hdr)

% Reads k- & kik-Net ascii files
% menandrin@gmail.com, 130528
%tic
if (nargin<3)
    o_hdr = 0;
end

 
fid    = fopen(traceFullName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

% End of line
eolLim = '\r?\n';

% Line
LinePattern = ['[^\r\n]*', eolLim];
thelines    = regexp(StrRay,LinePattern,'match');


% Read header lines
if (~o_hdr)
    meta = [];
else
    hdrLines         = thelines(1:nhdr);
    
    meta.eqTime      = hdrLines{1}(19:37);
    meta.eqLat       = str2num(hdrLines{2}(19:end));
    meta.eqLon       = str2num(hdrLines{3}(19:end));
    meta.eqZ         = str2num(hdrLines{4}(19:end));   % Eq Depth [km]
    meta.m           = str2num(hdrLines{5}(19:end));
    meta.sr          = str2num(hdrLines{11}(19:21));
    meta.stLat       = str2num(hdrLines{7}(19:end));
    meta.stLon       = str2num(hdrLines{8}(19:end));
    meta.stZ         = str2num(hdrLines{9}(19:end));   % Station height [m]
    meta.trec        = hdrLines{10}(19:end-1);   
    meta.scaleFactor = hdrLines{14}(19:end-1);
    meta.maxAcc      = str2num(hdrLines{15}(19:end));
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