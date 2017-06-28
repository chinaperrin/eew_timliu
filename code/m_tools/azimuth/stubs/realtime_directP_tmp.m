% OVERVIEW
% 1. Load test seismogram
% 2. Find p-wave onset on vertical comp
% 3. Cut out nsec seconds after p-pick
% 4. Run snippet through filter bank --> nband seismograms
% 5. Rotate seismograms about vertical axis, by <theta> degrees, and compute max-amps 
% 6. Find rotation that maximises high freqs on one and low freqs on the other axis;
%    the former corresponds to the radial, the latter to the transversal component
% 7. Measure and interpret direct P amplitude


% QUESTIONS
% - use velocity? displacement?
% - can i do filtering first and then rotate, or do i need to filter all rotated 
%   seismograms individually?


% TASKS
% - compare to other methods

%% 
clear all

addpath(genpath('../'))
addpath(genpath('../m_tools/'))
addpath(genpath('../m_picker/'))
addpath(genpath('../m_filter/'))
addpath(genpath('~/programs/matlab/fct_downloads/'))
addpath(genpath('~/programs/matlab/fct_downloads/subtightplot/'))

% Picker settings
fLow_px        = 3;
fUp_px         = 10; 
pxThreshold    = 5;
pxShift        = 15;             % Shift automatic px backward in time by <pxShift> samples 

o_plot_HvsV = 0;

% Filter specifications
fLow_prefilt = 1/20;

realisationType = 8;        % Choose filter implementation --> butter_pass_tdomain_f.m

snpLength = 0.5;
npad      = 100;
ntap      = 100;
T         = 40; 
tFrac     = 0.4;

vp             = 6;              % constant p-phase speed
vs             = 3.4;            % constant s-phase speed
sdelay         = (1/vs-1/vp);    % Delay of s-phase wrt/ p-phase per km [s/km]

ftSize = 12;

dataSetNames = {'/scratch/memeier/data/japan/k_kik/M6p_D25km/out/9bands_i27/'};
[TraceList,SkipList,fc,snpLength] = import_fbOut(dataSetNames);
nbands    = numel(fc(:,1));     

%%
idx_Z = find(strcmp(TraceList.orntCode,'Z'));
zList = TraceList.selectSubList(idx_Z);

IDX = find(zList.m>6 & zList.hypDist<40);
%idx = IDX(140);
ntr = numel(IDX);


for itr = 1:ntr
    
    idx = (IDX(itr));
    
    sr    = zList.sRate(idx);
    m     = zList.m(idx);
    hd    = zList.hypDist(idx);
    eqLat = zList.eqLat(idx);
    eqLon = zList.eqLon(idx);
    eqZ   = zList.eqZ(idx);
    stLat = zList.stationLat(idx);
    stLon = zList.stationLon(idx);
    
    %%
    
    % Azimuth of line from eq to station
    az = round(azimuth(eqLat,eqLon,stLat,stLon));
    
    %figure(131); clf; hold on
    %plot(eqLon,eqLat,'pk','markerFaceColor','y','markerSize',8)
    %plot(stLon,stLat,'vk','markerFaceColor','c','markerSize',8)
    
    %%
    zFullName           = zList.fullName{idx};
    [idx_Z,idx_E,idx_N] = find_corec_idx(zFullName,TraceList);
    zFullNameBak        = TraceList.fullName{idx_Z};
    eFullName           = TraceList.fullName{idx_E};
    nFullName           = TraceList.fullName{idx_N};
    
    {zFullName,zFullNameBak,eFullName,nFullName}'
    
    zTraceName = get_recordName(zFullName);
    fprintf(1,['\t',zTraceName,'\tsr=',num2str(sr),'Hz, m=',num2str(m),', @hD=', ...
        num2str(hd,'%4.0f'),'km, depth: ',num2str(eqZ),'km\n'])
    
    % Load example waveform
    [zraw,meta] = read_any_trace_noproc(zFullName);
    [eraw,~]    = read_any_trace_noproc(eFullName);
    [nraw,~]    = read_any_trace_noproc(nFullName);
    
    % zraw = zraw(1000:5000);
    % eraw = eraw(1000:5000);
    % nraw = nraw(1000:5000);
    
    ns   = numel(zraw);
    sr   = meta.sr;
    dt   = 1/sr;
    t    = dt:dt:ns*dt;
    
    % Remove early mean
    zm = zraw - mean(zraw(1:200));
    em = eraw - mean(eraw(1:200));
    nm = nraw - mean(nraw(1:200));
    
    ztp = taper(zm,sr,ntap);                                            % Taper waveform ...
    etp = taper(em,sr,ntap);                                            % Taper waveform ...
    ntp = taper(nm,sr,ntap);                                            % Taper waveform ...
    
    [z,~,~] = butter_pass_tdomain_f(ztp,fLow_prefilt,999,sr,7,npad,0);     % Prefilter ...
    [e,~,~] = butter_pass_tdomain_f(etp,fLow_prefilt,999,sr,7,npad,0);     % Prefilter ...
    [n,~,~] = butter_pass_tdomain_f(ntp,fLow_prefilt,999,sr,7,npad,0);     % Prefilter ...
    
    % Get pick
    [ppxIdx,snr] = sta_lta_picker_ntap(z,1/sr,pxThreshold,pxShift,1,ntap,1);
    1+1;
%     tppx   = t(ppxIdx);
%     dtsp   = hd*sdelay;         % ts - tp [s]
%     tspx   = tppx + dtsp;
%     spxIdx = ppxIdx + dtsp*sr;
%     
%     st   = tppx-1;        % Plot boundaries
%     et   = tspx+1;
%     sIdx = ppxIdx-1*sr;
%     eIdx = spxIdx+1*sr;
%     
%     gcf; 
%     
%     subplot(5,1,4:5); hold on
%     ymax = 1.1*max(abs(z));
%     line([t(spxIdx) t(spxIdx)],[-ymax ymax],'color',[0 .4 0],'lineWidth',2)
%     
%     % Run it through filter bank
%     [~,~,zout] = filterBank_1trace(z,sr,ppxIdx,snpLength,8,fc,npad);
%     [~,~,eout] = filterBank_1trace(e,sr,ppxIdx,snpLength,8,fc,npad);
%     [~,~,nout] = filterBank_1trace(n,sr,ppxIdx,snpLength,8,fc,npad);
%     
%     %plotTrace_minimal(z,zout,sr,fc)
%     %plotTrace_minimal(e,eout,sr,fc)
%     %plotTrace_minimal(n,nout,sr,fc)
%     
%     
%     
%     %hr   = abs(etmp)./abs(ztmp);
%     %hrc  = cumsum(hr,2);
%     
%     %%
%     
%     if (o_plot_HvsV)
%         figure(76); clf;
%         subplot(3,1,1); hold on
%         
%         ymax = 1.5*max(abs(e(sIdx:eIdx)));
%         set(gca,'xlim',[st et],'ylim',[-ymax ymax])
%         
%         
%         plot(t,z+ymax/2,'k')
%         plot(t,e)
%         plot(t,n-ymax/2,'r')
%         line([tppx tppx],[-ymax ymax],'color','r','linewidth',2)
%         line([tspx tspx],[-ymax ymax],'color',[0 0.4 0],'linewidth',2)
%         
%         band = 1;
%         
%         ztmp = zout{band};
%         etmp = eout{band};
%         
%         ztmp(1) = mean(ztmp(2:10));
%         etmp(1) = mean(etmp(2:10));
%         
%         ecum = cumsum(abs(etmp));
%         zcum = cumsum(abs(ztmp));
%         
%         ratio = ecum./zcum;
%         diff  = ecum-zcum;
%         
%         subplot(3,1,2); hold on
%         plot(t,ecum,'r')
%         plot(t,zcum,'k')
%         ymax = 1.1*max([max(ecum(sIdx:eIdx)), max(zcum(sIdx:eIdx))]);
%         legend('H','V')
%         plot(t(ppxIdx),0,'pk','markerFaceColor','r','markerSize',8)
%         line([tppx tppx],[-ymax ymax],'color','r','linewidth',2)
%         line([tspx tspx],[-ymax ymax],'color',[0 0.4 0],'linewidth',2)
%         set(gca,'xlim',[st et],'ylim',[0 ymax])
%         
%         subplot(3,1,3); hold on
%         [AX,H] = plotyy(t,ratio,t,diff);
%         ymax_r = 1.1*max(ratio(sIdx:eIdx));
%         ymin_r = 0.9*min(ratio(sIdx:eIdx));
%         ymax_d = 1.1*max(diff(sIdx:eIdx));
%         ymin_d = 0.9*min(diff(sIdx:eIdx));
%         set(AX(1),'ylim',[ymin_r ymax_r])
%         set(AX(2),'ylim',[ymin_d ymax_d])
%         set(AX,'xlim',[st et]);
%         legend('H / V','H - V')
%         plot(t(ppxIdx),0,'pk','markerFaceColor','r','markerSize',8)
%         line([tppx tppx],[-ymax ymax],'color','r','linewidth',2)
%         line([tspx tspx],[-ymax ymax],'color',[0 0.4 0],'linewidth',2)
%         
%         % %% Does it make a difference whether I a) filter the seismogram and then
%         % %  rotate it or b) rotate it and then filter?
%         % band  = 4;
%         % theta = 350;
%         % R     = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
%         %
%         % % a. filter then rotate
%         % e1      = eout{band};
%         % n1      = nout{band};
%         % en1     = [e1, n1];
%         % en1rot  = en1*R;
%         % e1rot   = en1rot(:,1);
%         % n1rot   = en1rot(:,2);
%         %
%         % % b. rotate then filter
%         % en2    = [e, n];
%         % en2rot = en2*R;
%         % [~,~,tmp1] = filterBank_1trace(en2rot(:,1),sr,ppxIdx,snpLength,8,fc,npad);
%         % [~,~,tmp2] = filterBank_1trace(en2rot(:,2),sr,ppxIdx,snpLength,8,fc,npad);
%         % e2rot  = tmp1{band};
%         % n2rot  = tmp2{band};
%         %
%         % figure(9); clf; hold on; plot(e2rot,'k','lineWidth',2); plot(e1rot,'r')
%         % % ==> Yes, they are the same, order does not play a role
%     end
%     
%     
%     %% Rotate seismograms and measure max amps between p- & s-phases in each freq band
%     
%     dtheta = 10;
%     theta  = dtheta:dtheta:360;
%     nrot   = numel(theta);
%     
%     d90  = -270:180:270;
%     az90 = az + d90;
%     na   = numel(az90);
%     
%     eMax = zeros(nbands,nrot);
%     nMax = zeros(nbands,nrot);
%     
%     figure(23); clf;
%     %ha = tight_subplot(9,1,[0 0.3],[0.05 0.05]);
%     ha = tight_subplot(10,1,[0.02 0],[0.08 0.05],[0.12 0.1]);
%     whitebg(23,'k')
%     
%     axes(ha(1));
%     plot(t,n,'w')
%     ymax = 1.1*max(abs(n(sIdx:eIdx)));
%     set(gca,'xlim',[st et],'ylim',[-ymax ymax])
%     line([tppx tppx],[-ymax ymax],'color',[0 .4 0],'lineStyle','--','lineWidth',2)
%     line([tspx tspx],[-ymax ymax],'color',[0 .4 0],'lineStyle','--','lineWidth',2)
%     set(gca,'XAxisLocation','top')
%     
%     title('Max. abs. p-phase amplitudes [m/s]','fontSize',ftSize)    
%     
%     for iband = 1:nbands
%         
%         % Write filtered wforms between p- & s-pick into ns-by-2 matrix H
%         ntmp = nout{iband}(ppxIdx:spxIdx);
%         etmp = eout{iband}(ppxIdx:spxIdx);
%         H = [ntmp etmp];
%         
%         for irot = 1:nrot
%             
%             R = [cosd(theta(irot)) -sind(theta(irot)); ...
%                 sind(theta(irot))  cosd(theta(irot))];
%             
%             Hrot  = H*R;
%             
%             nMax(iband,irot) = max(abs(Hrot(:,1)));
%             eMax(iband,irot) = max(abs(Hrot(:,2)));
%             
%             1+1;
%         end
%         
%         axes(ha(iband+1));
%         plot(theta,nMax(iband,:),'x')
%         ymax = 1.1*max(nMax(iband,:));
%         ymin = 0.9*min(nMax(iband,:));
%         %plot(theta,eMax(iband,:),'d')
%         hold on;
%         line([az90; az90],[repmat(ymin,1,na); repmat(ymax,1,na)],'color','r','lineStyle','--')
%         line([az az]    ,[ymin ymax],'color','w','lineWidth',2)
%         if (iband>4)
%             ylabel([num2str(fc(iband,1),'%3.1f'),' - ',num2str(fc(iband,2),'%3.1f')],'fontSize',ftSize)
%         else
%             ylabel([num2str(fc(iband,1),'%3.0f'),' - ',num2str(fc(iband,2),'%3.0f')],'fontSize',ftSize)
%         end
%         set(gca,'ylim',[ymin ymax],'xlim',[0 360])
%         
%         if (iband==nbands)
%             xlabel('Rotation angle [deg]','fontSize',ftSize)
%             %set(gca,'xTickLabel',{0:30:330})
%         else
%             set(gca,'xTickLabel',{})
%         end
%     end
    1+1;
end
