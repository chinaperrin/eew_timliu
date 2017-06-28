clear all

% Determine RT-backazimuth with 
% ND-method
% SH-method
% Lockman and Allen 2005, BSSA

% Compare to true backazimuth, as computed with hypocenter corrdinates


addpath(genpath('../'))
addpath(genpath('../m_tools/'))
addpath(genpath('../../../../../../gmpe/'))
addpath(genpath('../../../../../matlab/fct_downloads/'))

ftSize = 12; 

% Socal 90ies
dataSetNames = {'/scratch/memeier/data/japan/k_kik/M6p_D25km/out/9bands_i27/'}; 
[TraceList,SkipList,fc,snpLength] = import_fbOut(dataSetNames);

% Print/plot a summary of kept and skipped traces
[summary] = summary_fB(TraceList,SkipList,1);

vp             = 6;              % constant p-phase speed
vs             = 3.4;            % constant s-phase speed
sdelay         = (1/vs-1/vp);    % Delay of s-phase wrt/ p-phase per km [s/km]

% Extract records from within 100km
idx = find(cellfun(@(x) ~isempty(x), regexp(TraceList.fullName,'KGS0029703261731')));
%idx       = find(TraceList.hypDist<50);
TraceList = TraceList.selectSubList(idx);

% Station Lists
stFileName         = '../../../VS/data/stations/eewvs-stations_allCali.hinv.fmtd2';
[stationList_cali] = import_stationlist_fB(stFileName);
stFileName         = '../../../VS/data/stations/japan/sitepub_all_en.txt';
stationList_jp     = import_k_kik_net_list(stFileName);


%% Going through TraceList
n = numel(TraceList.m);
for i = 1:n
    
    
    % 1. Load wform and meta data
    % ---------------------------
    [status,result] = unix(['ls ',TraceList.fullName{i}]);
    
    if (status==0)
        
        % File found
        traceFullName = TraceList.fullName{i};
        slash_idx     = regexp(traceFullName,'/');
        traceName     = traceFullName(slash_idx(end)+1:end);
     
        % Find corecs
        [idx_Z,idx_E,idx_N] = find_corec_idx(traceFullName,TraceList);
        
        ornt   = TraceList.orntCode{i};
        instr  = TraceList.instrCode{i};
        m      = TraceList.m(i);
        hD     = TraceList.hypDist(i);
        pgv    = TraceList.pgv(i);
        stLat  = TraceList.stationLat(i);
        stLon  = TraceList.stationLon(i);
        eqLat  = TraceList.eqLat(i);
        eqLon  = TraceList.eqLon(i);
        eqDate = TraceList.eqDate{i};

        fprintf(1,[num2str(i),'. out of ',num2str(n),' waveform-files found\n'])
        fprintf(1,['   ',num2str(m),'M @ ',num2str(hD),'km\n'])
        titleStr = strcat([traceName,': M',num2str(m),' @ ',num2str(hD,'%4.0f'),'km, on ',eqDate]);
     
        % Read traces
        [S,sraw,meta] = read_any_trace_proc(TraceList.fullName{idx_Z},TraceList);
        z             = S.dsp;
        zv            = S.vel;
        ppxIdx        = meta.ppxIdx;         % sraw is as read from sac-file. NaN were removed (warning message in stdout),
        tz            = meta.t;             % and amplitudes have been devided by 100 to turn [cm/s] into [m/s] or [cm/s^2] into [m/s^2]
        sr            = meta.sr;
        
        [S,sraw,meta] = read_any_trace_proc(TraceList.fullName{idx_E},TraceList);
        e             = S.dsp;
        ev            = S.vel;
        te            = meta.t;
        
        [S,sraw,meta] = read_any_trace_proc(TraceList.fullName{idx_N},TraceList);
        n             = S.dsp;
        nv            = S.vel;
        tn            = meta.t;
        
        [z,n,e,t,ppxIdx] = harmonise_signal_times(z,n,e,tz,te,tn,ppxIdx);
        clear tz tn te
        
        [Z,~,~]  = butter_pass_tdomain_f(z,3,999,TraceList.sRate(i),7,100,1);
        [E,~,~]  = butter_pass_tdomain_f(e,3,999,TraceList.sRate(i),7,100,1);
        [N,~,~]  = butter_pass_tdomain_f(n,3,999,TraceList.sRate(i),7,100,1);
        [Zv,~,~] = butter_pass_tdomain_f(zv,3,999,TraceList.sRate(i),7,100,1);
        [Ev,~,~] = butter_pass_tdomain_f(ev,3,999,TraceList.sRate(i),7,100,1);
        [Nv,~,~] = butter_pass_tdomain_f(nv,3,999,TraceList.sRate(i),7,100,1);
        
        % Theoretical s-phase pick
        tppx   = t(ppxIdx);
        dtsp   = hD*sdelay;         % ts - tp [s]
        tspx   = tppx + dtsp;
        spxIdx = ppxIdx + round(dtsp*sr);
        
        % Plot boundaries
        st   = tppx-1;        
        et   = tspx+1;
        sIdx = ppxIdx-1*sr;
        eIdx = spxIdx+1*sr;
        
        
        % 2. Compute true backazimuth
        % ---------------------------
        az  = round(azimuth(eqLat,eqLon,stLat,stLon));
        bAz = round(azimuth(stLat,stLon,eqLat,eqLon));
        fprintf(1,'q: Is back-azimuth sth else than azimuth +/- 180deg?\n')
        %[bAz,ray] = get_backAzimuth(eqLat,eqLon,stLat,stLon);
       

        
        % 3. Compute real time backazimuth
        % --------------------------------

        % a. Nicolas Deichmann
        tw = tspx-tppx;   % Time over which angles are estimated
        [azimuth,incidence]   = azimuth_incidence_org(Ev(ppxIdx:spxIdx),Nv(ppxIdx:spxIdx),Zv(ppxIdx:spxIdx),1/sr,1);
        [azimuth2,incidence2] = azimuth_incidence(N,E,Z,t,1/sr,ppxIdx,tw,1);
        
        % b. Lockman and Allen, 2005, BSSA
        alpha  = 0.5;
        nz     = numel(Z);
        theta1 = zeros(nz,1);
        R_ze   = zeros(nz,1);
        R_zn   = zeros(nz,1);
        ns     = sr*5;
        
        for is = ppxIdx:ppxIdx+ns;
            R_ze(is)   = alpha*R_ze(is-1) + Z(is)*E(is);
            R_zn(is)   = alpha*R_zn(is-1) + Z(is)*N(is);
            theta1(is) = rad2deg(pi + atan(R_ze(is)/R_zn(is)));
        end
        rtBaz1 = mean(theta1(ppxIdx:ppxIdx+ns));


        
        % c. Seismic Handler

        
        % back-azimuth errors
        dbaz1 = bAz - rtBaz1; 
        
        % 4. Plot 

        figure(246); clf; hold on

        % Time window
        sIdx   = round(ppxIdx - 2*sr);      % t-index of window-start
        eIdx   = round(ppxIdx + 5*sr);      % t-index of window-end
        ymax   = 1.05*max([max(abs(Z(sIdx:eIdx))),max(abs(N(sIdx:eIdx))),max(abs(E(sIdx:eIdx)))]);
        
        
        % Vertical  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        subplot(4,3,1:2); hold on
        plot(t,Z,'k')
        line([t(ppxIdx) t(ppxIdx)],[-ymax ymax],'color','r')
        ylabel('UD [m]','fontSize',ftSize)
        axis([t(sIdx) t(eIdx) -ymax ymax])
        
        % North/South -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        subplot(4,3,4:5); hold on
        plot(t,N,'k')
        line([t(ppxIdx) t(ppxIdx)],[-ymax ymax],'color','r')
        ylabel('NS [m]','fontSize',ftSize)
        axis([t(sIdx) t(eIdx) -ymax ymax])
        
        % East/West  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        subplot(4,3,7:8); hold on
        plot(t,E,'k')
        line([t(ppxIdx) t(ppxIdx)],[-ymax ymax],'color','r')
        ylabel('EW [m]','fontSize',ftSize)
        axis([t(sIdx) t(eIdx) -ymax ymax])
        xlabel('Time [sec]','fontSize',ftSize)
        
        % Map  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        %subplot(4,3,[3,6])
        gcf; clf;
        station_map(stationList_cali)
        dlat = 0.5; dlon = 0.5;
        set(gca,'xlim',[stLon-dlon,stLon+dlon],'ylim',[stLat-dlat,stLat+dlat])
        hold on
        h1 = plot(stLon,stLat,'vk','markerFaceColor','y','markerSize',15);
        h2 = plot(eqLon,eqLat,'pk','markerFaceColor','r','markerSize',15);
        legend([h1,h2],'Station','Epicenter')
        title(['True back-azimuth: ',num2str(bAz,'%3.0f'),' [deg]'])
        %arrow([stLon,stLat],[eqLon,eqLat],'width',3,'faceColor','r')
        
        dlat = stLat - eqLat;
        dlon = stLon - eqLon;
        dl   = sqrt(dlat*dlat + dlon*dlon);
        
        % True back-azimuth
        arrowPt.lon = stLon + dl/1.2*sin(deg2rad(bAz));
        arrowPt.lat = stLat + dl/1.2*cos(deg2rad(bAz));
        arrow([stLon,stLat],[arrowPt.lon,arrowPt.lat],'width',3,'faceColor','r')
        
        % Lockman and Allen, 2005, BSSA
        arrowPt1.lon = stLon + dl/2*sin(deg2rad(rtBaz1));
        arrowPt1.lat = stLat + dl/2*cos(deg2rad(rtBaz1));
        arrow([stLon,stLat],[arrowPt1.lon,arrowPt1.lat],'width',2,'faceColor','b')
            
        1+1;
        
    end
end

        