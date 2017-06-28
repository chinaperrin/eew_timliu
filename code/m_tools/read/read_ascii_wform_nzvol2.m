function [s,meta] = read_ascii_wform_nzvol2(traceFullName)

% Reads PROCESSED ascii waveform file from the strong motion data base in
% New Zealand. Meta info is read only for first block and assumed to be the
% same for other two blocks. This "includes the start of buffer time" which
% I think gives the absolute time of the first waveform sample.
%
% menandrin@gmail.com, 161118
 
nhdr = 26;
ncol = 10;

fid    = fopen(traceFullName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

eolLim      = '\r?\n';
LinePattern = ['[^\r\n]*', eolLim];
thelines    = regexp(StrRay,LinePattern,'match');



%% META INFO
ln2 = regexp(thelines{2},'\s+','split');
meta.station.name = ln2{2};
meta.station.lat  = -(str2double(ln2{3}) + str2double(ln2{4})/60 + str2double(ln2{5}(1:2))/3600);
meta.station.lon  =   str2double(ln2{6}) + str2double(ln2{7})/60 + str2double(ln2{8}(1:2))/3600;

ln5 = regexp(thelines{5},'\s+','split');
meta.instrument.period = str2double(ln5{5});
meta.instrument.damp   = str2double(ln5{8});

ln9 = regexp(thelines{9},'\s+','split');
meta.eq.lat = -(str2double(ln9{2}) + str2double(ln9{3})/60 + str2double(ln9{4}(1:2))/3600);
meta.eq.lon =   str2double(ln9{5}) + str2double(ln9{6})/60 + str2double(ln9{7}(1:2))/3600;
meta.eq.bearing = ln9{9};
meta.eq.epidist = str2double(ln9{11}(1:end-2));
meta.eq.depth   = str2double(ln9{13}(1:end-2));
meta.eq.Mw      = str2double(ln9{15});

ln10 = regexp(thelines{10},'\s+','split');
meta.record.ns       = str2double(ln10{4});
meta.record.duration = str2double(ln10{6});

ln11 = regexp(thelines{11},'\s+','split');
meta.record.sr     = 1/str2double(ln11{5});
meta.record.filter = thelines{12};

ln13 = regexp(thelines{13},'\s+','split');
meta.instrument.comp1 = ln13{2};

ln14 = regexp(thelines{14},'\s+','split');
ln15 = regexp(thelines{15},'\s+','split');
ln16 = regexp(thelines{16},'\s+','split');
meta.record.pga  = str2double(ln14{3});
meta.record.tpga = str2double(ln14{6});
meta.record.pgv  = str2double(ln15{3});
meta.record.tpgv = str2double(ln15{6});
meta.record.pgd  = str2double(ln16{3});
meta.record.tpgd = str2double(ln16{6});

% Origin time 
ln17 = regexp(thelines{17},'\s+','split');
ln18 = regexp(thelines{18},'\s+','split');
ln20 = regexp(thelines{20},'\s+','split');
yr   = str2double(ln17{2});
mt   = str2double(ln17{3});
dy   = str2double(ln17{4});
hr   = str2double(ln17{5});
mn   = str2double(ln17{6});
sc   = str2double(ln17{7});
t0   = datetime(yr,mt,dy,hr,mn,sc*1e-1,'InputFormat','yy,mm,dd,HH,MM,SS');
meta.eq.t0 = t0;


% Buffer start time
yr   = str2double(ln17{10});
mt   = str2double(ln17{11});
dy   = str2double(ln18{10});
hr   = str2double(ln18{11});
mn   = str2double(ln20{10});
sc   = str2double(ln20{11});
ts   = datetime(yr,mt,dy,hr,mn,sc*1e-3,'InputFormat','yy,mm,dd,HH,MM,SS');
meta.record.ts = ts;

% Samples concatenated to waveform
ns_prepended    = str2double(ln20{3});
dt_prepended    = ns_prepended/meta.record.sr;
dts             = etime(datevec(ts),datevec(t0))-dt_prepended;
meta.record.dts = dts;


%%  BLOCK 1 ________________________________________________________________
nsLine = regexp(thelines{20},'\s+','split');
na     = str2double(nsLine(5));
nv     = str2double(nsLine(6));
nd     = str2double(nsLine(7));
ns1    = nhdr+1;                  % First data line of block 1 
ne1    = nhdr+ceil(na/10)+ceil(nv/10)+ceil(nd/10);  % Last  data line of block 1

% Acceleration
nla    = floor(na/10);
nae    = ns1+nla-1;    % Last line of acc. wform
sraw   = str2num(cell2mat(thelines(ns1:nae)));
acc.h1 = reshape(sraw',nla*ncol,1);
nspare = rem(na,10);   % If last line is not fully populated ...
if nspare~=0; nae      = nae+1; 
              lastline = str2num(cell2mat(thelines(nae)))';
              acc.h1   = [acc.h1; lastline];
end
%clf; plot(acc.h1)

% Velocity
nlv    = floor(nv/10);
nvs    = nae+1;
nve    = nvs+nlv-1;    % Last line of acc. wform
sraw   = str2num(cell2mat(thelines(nvs:nve)));
vel.h1 = reshape(sraw',nlv*ncol,1);
nspare = rem(nv,10);   % If last line is not fully populated ...
if nspare~=0; nve      = nve+1; 
              lastline = str2num(cell2mat(thelines(nve)))';
              vel.h1   = [vel.h1; lastline];
end

% Displacement
nld    = floor(nd/10);
nds    = nve+1;
nde    = nds+nld-1;    % Last line of acc. wform
sraw   = str2num(cell2mat(thelines(nds:nde)));
dsp.h1 = reshape(sraw',nld*ncol,1);
nspare = rem(nd,10);   % If last line is not fully populated ...
if nspare~=0; nde      = nde+1; 
              lastline = str2num(cell2mat(thelines(nde)))';
              dsp.h1   = [dsp.h1; lastline];
end




%%  BLOCK 2 ________________________________________________________________
n02          = ne1+1;
isRightLine2 = ~isempty(regexp(thelines{n02},'GNS Science','start'));
ns2          = n02+nhdr;
ns2Line      = regexp(thelines{n02+19},'\s+','split');
na2          = str2double(ns2Line(5));
nv2          = str2double(ns2Line(6));
nd2          = str2double(ns2Line(7));
ne2          = ns2-1+ceil(na2/10)+ceil(nv2/10)+ceil(nd2/10);  % Last  data line of block 2

ln13 = regexp(thelines{n02+12},'\s+','split');
meta.instrument.comp2 = ln13{2};

% Acceleration
nla2   = floor(na2/10);
nae2   = ns2+nla2-1;    % Last line of acc. wform
sraw   = str2num(cell2mat(thelines(ns2:nae2)));
acc.h2 = reshape(sraw',nla2*ncol,1);
nspare = rem(na2,10);   % If last line is not fully populated ...
if nspare~=0; nae2      = nae2+1; 
              lastline = str2num(cell2mat(thelines(nae2)))';
              acc.h2   = [acc.h2; lastline];
end

% Velocity
nlv2   = floor(nv2/10);
nvs2   = nae2+1;         % First line of vel wform
nve2   = nvs2+nlv2-1;    % Last  line of vel wform
sraw   = str2num(cell2mat(thelines(nvs2:nve2)));
vel.h2 = reshape(sraw',nlv2*ncol,1);
nspare = rem(nv2,10);   % If last line is not fully populated ...
if nspare~=0; nve2      = nve2+1; 
              lastline = str2num(cell2mat(thelines(nve2)))';
              vel.h2   = [vel.h2; lastline];
end

% Displacement
nld2   = floor(nd2/10);
nds2   = nve2+1;         % First line of dsp wform
nde2   = nds2+nld2-1;    % Last  line of dsp wform
sraw   = str2num(cell2mat(thelines(nds2:nde2)));
dsp.h2 = reshape(sraw',nld2*ncol,1);
nspare = rem(nd2,10);   % If last line is not fully populated ...
if nspare~=0; nde2     = nde2+1; 
              lastline = str2num(cell2mat(thelines(nde2)))';
              dsp.h2   = [dsp.h2; lastline];
end



%%  BLOCK 3
n03          = ne2+1;
isRightLine3 = ~isempty(regexp(thelines{n03},'GNS Science','start'));
ns3          = n03+nhdr;
ns3Line      = regexp(thelines{n03+19},'\s+','split');
na3          = str2double(ns3Line(5));
nv3          = str2double(ns3Line(6));
nd3          = str2double(ns3Line(7));
%ne3          = ns3-1+na3/10+nv3/10+nd3/10;  % Last  data line of block 3

ln13 = regexp(thelines{n03+12},'\s+','split');
meta.instrument.comp3 = ln13{2};

% Acceleration
nla3   = floor(na3/10);
nae3   = ns3+nla3-1;    % Last line of acc. wform
sraw   = str2num(cell2mat(thelines(ns3:nae3)));
acc.z  = reshape(sraw',nla3*ncol,1);
nspare = rem(na3,10);   % If last line is not fully populated ...
if nspare~=0; nae3      = nae3+1; 
              lastline = str2num(cell2mat(thelines(nae3)))';
              acc.z   = [acc.z; lastline];
end

% Velocity
nlv3   = floor(nv3/10);
nvs3   = nae3+1;         % First line of vel wform
nve3   = nvs3+nlv3-1;    % Last  line of vel wform
sraw   = str2num(cell2mat(thelines(nvs3:nve3)));
vel.z  = reshape(sraw',nlv3*ncol,1);
nspare = rem(nv3,10);   % If last line is not fully populated ...
if nspare~=0; nve3     = nve3+1; 
              lastline = str2num(cell2mat(thelines(nve3)))';
              vel.z   = [vel.z; lastline];
end

% Displacement
nld3   = floor(nd3/10);
nds3   = nve3+1;         % First line of dsp wform
nde3   = nds3+nld3-1;    % Last  line of dsp wform
sraw   = str2num(cell2mat(thelines(nds3:nde3)));
dsp.z  = reshape(sraw',nld3*ncol,1);
nspare = rem(nd3,10);   % If last line is not fully populated ...
if nspare~=0; nde3      = nde3+1; 
              lastline = str2num(cell2mat(thelines(nde3)))';
              dsp.z   = [dsp.z; lastline];
end




%% 
if ~isRightLine2 ||~isRightLine3; error('Line numbers screwed up somehow. check. now.'); end

s.acc = acc;
s.vel = vel;
s.dsp = dsp;