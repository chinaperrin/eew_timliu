%% Read xls file that contains finite fault specification; 
%  File has been provided and updated by ABaltay, 160904
%
%  Some large events for which there are waveforms in the traceList data
%  base are missing in this finite fault source data base:
%  . 2016 Kumamoto, 20160416012500
%  . 20110411171600, M6.7, Tohoku aftershock
%  . 20000730212500, M6.5, Izu Islands

clear all
addpath(genpath('../'))
addpath(genpath('../../../matlab/'))

load /scratch/memeier/fbout/i38/causal_0p5/gmpp_targetM6p_trainM4p/out_tr64814_pid6802_all2242.mat
zList      = out.zList;
eqs        = zList.prop.eqs;
ntr        = numel(zList.m);
clear out

o.printPng    = 1;
o.saveOut     = 1;
o.outFullName = '~/Documents/eq/data/faults/new/allFaults_i38';

%% Read Annemarie Baltays' excel-file 
xlsFullName = '~/Documents/eq/data/faults/NGA_crossCheck_ii.xlsx';
[num1, str1]=xlsread(xlsFullName,'Sheet2', 'A1:Z30');
[num2, str2]=xlsread(xlsFullName,'Sheet3');


nlines = 25;
neqs   = nlines-1;

% Double entries:
% Landers:      LANDERS_1992 / 3031111  (entires in xls files with eventID as eqname were manually deleted)
% Hector Mine:  HECTOR_1999  / 9108652
% Northridge:   NORTHR_1994  / 3144585

% "A multi-segment rupture consists of two or more contiguous planar quadrilaterals 
% joined along their down-dip edges, sharing a single hypocenter."
% "A multi-fault rupture occurs on two or more non-contiguous surfaces (which might 
% each be multi-segment) and each surface having its own hypocenter."
% ==> I concatenate all segments, irrespective of whether they are on same fault or not
kocaeli.faultlines    = 1;
landers.faultlines    = 2;
duzce.faultlines      = 3;
chichi.faultlines     = 4;
denali.faultlines     = 6:8;            % --> Line numbers on Sheet3
kobe.faultlines       = 10:11;
elmayor.faultlines    = 13:16;
wenchuan.faultlines   = 18:19;
hectormine.faultlines = 21:23;

isSimpleModel                     = true(neqs,1);
isSimpleModel([1:6,8,9,10,11,12]) = 0;

nFaults          = ones(neqs,1);
nFaults([19,23]) = 0;
nFaults(1)  = numel(denali.faultlines);
nFaults(2)  = numel(wenchuan.faultlines);
nFaults(3)  = numel(chichi.faultlines);
nFaults(4)  = numel(kocaeli.faultlines);
nFaults(5)  = numel(landers.faultlines);
nFaults(6)  = numel(landers.faultlines);
nFaults(8)  = numel(elmayor.faultlines);
nFaults(9)  = numel(hectormine.faultlines);
nFaults(10) = numel(hectormine.faultlines);
nFaults(11) = numel(duzce.faultlines);
nFaults(12) = numel(kobe.faultlines);

nameCol   = 7;
strikeCol = 11;
dipCol    = 12;
lengthCol = 13;
widthCol  = 14;
rakCol    = 15;
clatCol   = 16;     % Upper left corner
clonCol   = 17;
czCol     = 18;
mCol      = 2;
latCol    = 3;      % Hypocenter (?)
lonCol    = 4;
zCol      = 5;

fsrc.eqname   = cell(nlines-1,1);
fsrc.m        = zeros(nlines-1,1);
fsrc.lat      = zeros(nlines-1,1);
fsrc.lon      = zeros(nlines-1,1);
fsrc.z        = zeros(nlines-1,1);
fsrc.segments = cell(nlines-1,1);

for iline = 2:nlines
    ieq    = iline-1;
    eqName = str1{iline,nameCol};
    if isempty(eqName); eqName = num2str(num1(iline-1,nameCol)); end
    eqName           = strrep(eqName,' ','');
    fsrc.eqname{ieq} = eqName;

    isDenali     = strcmp(eqName,'DENALI_2002');
    isWenchuan   = strcmp(eqName,'wenchuan');
    isChichi     = strcmp(eqName,'CHICHI00_1999');
    isKocaeli    = strcmp(eqName,'KOCAELI_1999');
    isLanders    = strcmp(eqName,'LANDERS_1992');
    isDuzce      = strcmp(eqName,'DUZCE_1999');
    isHectormine = strcmp(eqName,'HECTOR_1999');
    isElmayor    = strcmp(eqName,'14607652');
    isKobe       = strcmp(eqName,'KOBE_1995');
    isKumamoto   = strcmp(eqName,'20160416012500');
    
    if strcmp(eqName,'3031111'); isLanders=1;    iline=iline+1; end
    if strcmp(eqName,'9108652'); isHectormine=1; iline=iline-1; end
    if strcmp(eqName,'3144585'); isNorthRidge=1; iline=iline-1; end
    
    fsrc.m     (ieq) = num1(iline-1,mCol);
    fsrc.lat   (ieq) = num1(iline-1,latCol);
    fsrc.lon   (ieq) = num1(iline-1,lonCol);
    fsrc.z     (ieq) = num1(iline-1,zCol);
    
    % Concatenate all segments, irrespective of whether they are on same
    % fault or not
    segments.clat   = zeros(100,1);
    segments.clon   = zeros(100,1);
    segments.cz     = zeros(100,1);
    segments.strike = zeros(100,1);
    segments.dip    = zeros(100,1);
    segments.length = zeros(100,1);
    segments.width  = zeros(100,1);
    ctSegment       = 0;
    
    for ifault = 1:nFaults(ieq)

        % Read line: if is simple model, read line from Sheet2, else read
        % one of the nFaults(ieq) lines in Sheet3
        if isSimpleModel(ieq); 
           theline = [99999,99999,num1(iline-1,:)]; 
            nseg   = 1;
        else
            if     isDenali;     lnNumber = denali.faultlines(ifault);
            elseif isChichi;     lnNumber = chichi.faultlines(ifault);
            elseif isWenchuan;   lnNumber = wenchuan.faultlines(ifault);
            elseif isKocaeli;    lnNumber = kocaeli.faultlines(ifault);
            elseif isLanders;    lnNumber = landers.faultlines(ifault);
            elseif isDuzce;      lnNumber = duzce.faultlines(ifault);
            elseif isHectormine; lnNumber = hectormine.faultlines(ifault);
            elseif isElmayor;    lnNumber = elmayor.faultlines(ifault);
            elseif isKobe;       lnNumber = kobe.faultlines(ifault);
            end
            theline = num2(lnNumber,:);
            nseg    = theline(12);
        end
        
        
        for iseg = 1:nseg
            
            ctSegment = ctSegment+1;
            
            % Read segment info
            segments.clat  (ctSegment) = theline(19+(iseg-1)*9);
            segments.clon  (ctSegment) = theline(20+(iseg-1)*9);
            segments.cz    (ctSegment) = theline(21+(iseg-1)*9);
            segments.strike(ctSegment) = theline(14+(iseg-1)*9);
            segments.dip   (ctSegment) = theline(15+(iseg-1)*9);
            segments.length(ctSegment) = theline(16+(iseg-1)*9);
            segments.width (ctSegment) = theline(17+(iseg-1)*9);
        end
        1+1;
    end
    
    % Truncate matrices
    segments.clat   = segments.clat  (1:ctSegment);
    segments.clon   = segments.clon  (1:ctSegment);
    segments.cz     = segments.cz    (1:ctSegment);
    segments.strike = segments.strike(1:ctSegment);
    segments.dip    = segments.dip   (1:ctSegment);
    segments.length = segments.length(1:ctSegment);
    segments.width  = segments.width (1:ctSegment);
    
    if isKumamoto
        
        % This is a file that Han Yue has sent me on 160904; it contains
        % point coordinates of surface rupture observations, i.e. these are
        % not fault segment corners. But as long as I assume short
        % segments, this should only result in small fault distance errors.
        kumamotoData = importdata('~/Documents/eq/data/faults/2016_kumamoto.txt');
        kumamotoData = kumamotoData(~isnan(kumamotoData(:,1)),:);
        
        segments.clat   = kumamotoData(:,1);
        segments.clon   = kumamotoData(:,2);
        segments.cz     = kumamotoData(:,3);
        segments.strike = 224*ones(size(segments.clat));   % from http://earthquake.usgs.gov/earthquakes/eventpage/us20005iis#finite-fault
        segments.dip    = 90*ones(size(segments.clat));      
        segments.length = 1*ones(size(segments.clat));     % from trial and error
        segments.width  = 15*ones(size(segments.clat));    % from trial and error
    end
    
    fsrc.segments{ieq} = segments;
end

fsrc.z(fsrc.z<0)=0;

hd_bak = zList.hypDist;
allFaults.eqname   = cell(neqs,1);
allFaults.segments = cell(neqs,1);
for ieq = 1:neqs

    fig.plot = 0;
    fig.view = [20 80];
    
    % Construct point cloud fault model for each segment of fault, then
    % concatenate
    nseg      = numel(fsrc.segments{ieq}.clat);
    fault.lat = [];
    fault.lon = [];
    fault.x   = [];
    fault.y   = [];
    fault.z   = [];
    for iseg = 1:nseg
        thisSegment = fsrc.segments{ieq}; 
        clat   = thisSegment.clat(iseg); 
        clon   = thisSegment.clon(iseg); 
        cz     = thisSegment.cz(iseg);
        strike = thisSegment.strike(iseg);
        dip    = thisSegment.dip(iseg);
        L      = 0:1:thisSegment.length(iseg);
        W      = (0:1:thisSegment.width(iseg))';
        newfault  = get_rectangular_fault_corners(clat,clon,cz,strike,dip,L,W,fig);
        
        fault.lat = [fault.lat; newfault.lat(:)];
        fault.lon = [fault.lon; newfault.lon(:)];
        fault.x   = [fault.x  ; newfault.x(:)];
        fault.y   = [fault.y  ; newfault.y(:)];
        fault.z   = [fault.z  ; newfault.z(:)];
    end
    
    if isnan(clat) |isnan(clon) |isnan(cz); % There is one entry where corner coordinates are missing
        fault.lat=[]; fault.lon=[]; fault.x=[]; fault.y=[]; fault.z=[];
    end


    allFaults.eqname{ieq}   = fsrc.eqname{ieq};
    allFaults.segments{ieq} = fault;
    
    ii    = find(strcmp(fsrc.eqname{ieq},eqs.name) & abs(fsrc.m(ieq)-eqs.m)<.3);
    trIdx = eqs.traceId{ii};
    stLat = zList.stationLat(trIdx);
    stLon = zList.stationLon(trIdx);
    stAlt = zList.stationAlt(trIdx);
    
    figure(812); clf; hold on; box on; grid on;
    plot3(fault.lon,fault.lat,fault.z,'xk','markerSize',5)
    set(gca,'zdir','reverse')
    plot3(stLon,stLat,zeros(size(stLon)),'vk','markerFaceColor','y','markerSize',15)
    view([10 60])
    xlm = get(gca,'xlim');
    ylm = get(gca,'ylim');
    zlm = get(gca,'zlim');
    fill3([xlm(1) xlm(2) xlm(2) xlm(1)],[ylm(1) ylm(1) ylm(2) ylm(2)],[zlm(1) zlm(1) zlm(1) zlm(1)],'k','faceAlpha',.1)
    title(strrep(eqs.name{ii},'_',' '))
    
    if o.printPng
        set(gcf,'PaperPositionMode','auto');
        pngFullName = sprintf('~/programs/filterBank/fig/i38/finSrc/new/%s_m%3.1f',fsrc.eqname{ieq},fsrc.m(ieq));
        pngFullName = strrep(pngFullName,'.','p');
        print('-dpng',pngFullName)
    end
1+1;
end


if o.saveOut
    save(o.outFullName,'allFaults')
end