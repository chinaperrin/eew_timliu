function [respList_SF] = import_SAFOD_simpleResp_list(listFullName)
% listFullName = '/scratch/memeier/data/nocal/safod/stations/SF_simple_responses.txt'

fid    = fopen(listFullName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

% End of line
eolLim = '\r?\n';

% Line
LinePattern = ['[^\r\n]*', eolLim];
thelines    = regexp(StrRay,LinePattern,'match');
bodyLines   = thelines(2:end);
nlines      = numel(bodyLines);
%headerLine = thelines(1);

% Initiate arrays
net     = cell (nlines,1);
sta     = cell (nlines,1);
cha     = cell (nlines,1);
loc     = cell (nlines,1);
ondate  = cell (nlines,1);
offdate = cell (nlines,1);
natFreq = zeros(nlines,1);
hpCfreq = zeros(nlines,1);
lpCfreq = zeros(nlines,1);
damping = zeros(nlines,1);
gain    = zeros(nlines,1);
azi     = zeros(nlines,1);
dip     = zeros(nlines,1);
unit    = cell (nlines,1);

for iline=1:nlines
    
    ctLine = bodyLines{iline};
    fields = regexp(ctLine,'\s+', 'split');
    
    net(iline)     = fields(1);
    sta(iline)     = fields(2);
    cha(iline)     = fields(3);
    loc(iline)     = fields(4);
    ondate(iline)  = fields(5);
    offdate(iline) = fields(6);
    natFreq(iline) = str2double(fields{7});
    hpCfreq(iline) = str2double(fields{8});
    lpCfreq(iline) = str2double(fields{9});
    damping(iline) = str2double(fields{10});
    gain(iline)    = str2double(fields{11});
    azi(iline)     = str2double(fields{12});
    dip(iline)     = str2double(fields{13});
    unit(iline)    = fields(14);
end

% Write to table
respList_SF         = table;
respList_SF.net     = net;
respList_SF.sta     = sta;
respList_SF.cha     = cha;
respList_SF.loc     = loc;
respList_SF.ondate  = ondate;
respList_SF.offdate = offdate;
respList_SF.natFreq = natFreq;
respList_SF.hpCfreq = hpCfreq;
respList_SF.lpCfreq = lpCfreq;
respList_SF.damping = damping;
respList_SF.gain    = gain;
respList_SF.azi     = azi;
respList_SF.dip     = dip;
respList_SF.unit    = unit;

%respList_SF(end-10:end,:)