function pickList = read_shakeAlert_pick_log_files(pickFileFullName)

% Reads pick files that JA generates when downloading waveforms which have
% an STA/LTA trigger in the ShakeAlert logfiles.
%
% mam, 170411


fid    = fopen(pickFileFullName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

eolLim      = '\r?\n';
LinePattern = ['[^\r\n]*', eolLim];
thelines    = regexp(StrRay,LinePattern,'match');
nlines      = numel(thelines);

pickList.network        = cell(nlines,1);
pickList.station        = cell(nlines,1);
pickList.channel        = cell(nlines,1);
pickList.location       = cell(nlines,1);
pickList.datetimestring = cell(nlines,1);
pickList.datetime       = cell(nlines,1);

for iline=1:nlines

    print_iteration_numbers(iline,nlines,'thousands')
    
    thisline = thelines{iline};
    thisline = strrep(thisline,'Peak',' Peak');
    fields   = regexp(thisline,' ','split');

    pickList.network{iline}  = fields{2};
    pickList.station{iline}  = fields{3};
    pickList.channel{iline}  = fields{4};
    pickList.location{iline} = strrep(fields{5},'''','');
    
    dtstring                       = fields{14};
    pickList.datetimestring{iline} = dtstring(1:end-1);
    pickList.datetime{iline}       = datetime(pickList.datetimestring{iline},'InputFormat','yyyy/MM/dd,HH:mm:ss.SSSS');
end