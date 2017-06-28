function [sraw,meta] = read_ascii_wform_tk(traceFullName,o_hdr,o_wf)

%traceFullName = '/scratch/memeier/VS/data/wform/turkey/19990817/19990817000151_5401.txt';
%nhdr = 18;

% Reads k- & kik-Net ascii files
% menandrin@gmail.com, 130528
%tic
if (nargin<2)
    o_hdr = 0;
    o_wf  = 1;
end

nhdr = 18;

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
    
    hdrLines      = thelines(1:nhdr);
    
    meta.eqDate   = regexp(hdrLines{3},'[0-9]+/[0-9]+/[0-9]+','match');
    meta.eqTime   = regexp(hdrLines{3},'[0-9]+:[0-9]+:([0-9]*[.,]?[0-9]+)','match');
    coords        = regexp(hdrLines{4},'\d*\.\d*','match');
    meta.eqLat    = str2double(coords{1});
    meta.eqLon    = str2double(coords{2});
    meta.eqZ      = str2double(cell2mat(regexp(hdrLines{5},'[-+]?([0-9]*[.,]?[0-9]+)','match')));
    meta.m        = str2double(cell2mat(regexp(hdrLines{6},'([0-9]*[.,]?[0-9]+)','match')));
    
    sInt          = str2double(cell2mat(regexp(hdrLines{14},'([0-9]*[.,]?[0-9]+)','match')));
    meta.sr       = 1/sInt;
    
    coords        = regexp(hdrLines{8},'([0-9]*[.,]?[0-9]+)','match');
    meta.stLat    = str2double(coords{1});
    meta.stLon    = str2double(coords{2});
    meta.stAlt    = str2double(cell2mat(regexp(hdrLines{9},'[0-9]*[,.]?[0-9]+','match')));   % Station height [m]
    if (isnan(meta.stAlt)); meta.stAlt = 0; end
    meta.stID     = regexp(hdrLines{7},'\d*','match');
    
    meta.trec     = regexp(hdrLines{12},'\d*:\d*:\d*','match');   
    meta.trecDate = regexp(hdrLines{12},'\d*/\d*/\d\d\d\d','match');   
    
    meta.pga      = regexp(hdrLines{15},'([0-9]*[.,]?[0-9]+)','match');
end



% Read header lines
if (~o_wf)
    sraw = [];
else
    % Read waveform amplitudes
    thelines = thelines(nhdr+1:end);
    nlines   = numel(thelines);
    raw      = str2num(cell2mat(thelines));
    sraw.ns  = raw(:,1);
    sraw.ew  = raw(:,2);
    sraw.ud  = raw(:,3);
end