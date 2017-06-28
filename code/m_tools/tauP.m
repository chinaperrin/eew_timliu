function [out] = tauP(traceFullName,trList,tauP0,alpha,o_plot,ts,te)

% Computes tauP as defined in Allen and Kanamori 2003. 
% 
% Olson and Allen, 2005    use lowpass filter with 3Hz, alpha=.99, t=0.05-4sec
% Allen and Kanamori, 2004 use lowpass fiter with 3 or 10Hz 
% menandrin@gmail.com, 130730

% Some parameters that different authors tweak in different ways
% - corner frequency of lowpass filter  <fUp>
% - time window for tau_p estimation    <tWindow>
% - alpha smoothing factor              <alpha>

% INPUT
% - tauP0       tauP estimation window ends at <tauP0> seconds after p-pick
% - alpha       smoothing factor

% OUTPUT
% tauP_out.m_l        : m-estimation with low-m  formula from Allen and Kanamori 2003  
% tauP_out.m_h        : m-estimation with high-m formula from Allen and Kanamori 2003
% tauP_out.m_olson    : m-estimation with formula from Olson and Allen, 2005
% tauP_out.tauP_max_l : tauP-max value from Allen and Kanamori 2003, T = 2sec, lowpass with fc = 10Hz;
% tauP_out.tauP_max_h = tauP-max value from Allen and Kanamori 2003, T = 4sec, lowpass with fc = 3Hz;

% NOTES
% - 8ung for NF records: if the S-wave is contained in the evaluation
%   window, the magnitude is typically strongly overestimated
% - The tauP formula is independent of the unit of the waveforms
% - Notes from R. Allen's article in EEW book: "for 100 sps data alpha=0.99,
%   for 20 sps data alpha=0.95"
% -------------------------------------------------------------------------


global vs

%% 1. Preparations

% Lowpass filter
fUp    = 3;         % Corner frequency of lowpass filter
fOrder = 4;
fMode  = 'causal'; 

startTime = 0.5;  % tauP estimation window starts <startTime> seconds after p-pick


if (o_plot)
    fprintf(1,sprintf('\n   Note: tauP is computed as the maximum between %3.1f - %3.1f sec after the p-pick\n',startTime,tauP0))
end

% Plot settings
lnWidth = 1;
ftSize  = 15;
if (nargin<4); o_plot = false; end

% Find entry in traceList  
idx = find(strcmp(traceFullName,trList.fullName));

% Read meta-data
m         = trList.m(idx);
r         = trList.hypDist(idx);
orntCode  = trList.orntCode{idx};

% Avoid using the S-phase: 
if r<vs*tauP0 
    tauP0 = r/vs;
    fprintf(1,'tauPmax window truncated to avoid using S-phase\n')
end

if (~strcmp(orntCode,'Z'))
    fprintf(1,'Warning: not a vertical component! This is not how tau_p is defined ...\n')
end

% Read and pre-process waveforms (Acc/Vel/Dsp)
[S,meta] = read_any_trace(trList.fullName{idx},trList,1);
acc      = S.acc;
vel      = S.vel;
t        = meta.t;
ppxIdx   = meta.ppxIdx;
sr       = meta.sr;
tppx     = t(ppxIdx);

% Select time interval
% if (nargin>4)
%     sIdx  = ppxIdx - tmin*sr;
%     eIdx  = ppxIdx + tmax*sr;
%     acc   = acc(sIdx:eIdx);
%     vel   = vel(sIdx:eIdx);
%     t     = t  (sIdx:eIdx);
%     ppxIdx = ppxIdx - (sIdx-1);
% end




%% 2. m_l
%  Allen and Kanamori 2003 use this estimate for events with 3.0<m<5.0
%  Filter broadband data with 10Hz lowpass filter. Use 2s (or 1s) of wform.
%  m_l = 6.3*log(Tp_max)+7.1

% Apply lowpass filter to both time series
[acc] = bworth(acc,sr,fUp,'low',fOrder,fMode);
[vel] = bworth(vel,sr,fUp,'low',fOrder,fMode);
n     = length(vel);

% Compute tau_p from p-pick on
tau_p = zeros(1,n);
Vi    = zeros(1,n);
Ai    = zeros(1,n);

%for i = ppxIdx+1:n
for i = 2:n
    Vi(i)    = alpha*Vi(i-1) + vel(i)^2;
    Ai(i)    = alpha*Ai(i-1) + acc(i)^2;
    tau_p(i) = 2*pi*sqrt(Vi(i)/Ai(i));
end

% Find max-value of tauP_l in first 0.5 - 4 seconds since p-pick
sIdx          = ppxIdx + startTime*sr;
eIdx          = ppxIdx + round(tauP0*sr)-1;
if eIdx>n; eIdx=n; end
if sIdx<1; sIdx=1; end
[tauPmax,idx] = max(tau_p(sIdx:eIdx));
idxMax        = idx+sIdx;               % Index at which tauPmax is measured

% Compute magnitudes with regression coeffs from literature
out.m_allen = 6.3 *log10(tauPmax) + 7.1;
out.m_olson = 7.14*log10(tauPmax) + 5.93;
out.tauPmax = tauPmax;
out.idxMax  = idxMax;


% %% 3. m_h
% %  Allen and Kanamori 2003 use this estimate for events with 4.5<m
% %  Filter broadband data with 3Hz lowpass filter, using 4s
% %  m_h = 7.0*log(Tp_max)+5.9
% 
% % Apply lowpass filter to both time series
% fUp         = 3; 
% [acc_h,~,~] = butter_pass_tdomain_f(acc,999,fUp,sr,6,npad,0);
% [vel_h,~,~] = butter_pass_tdomain_f(vel,999,fUp,sr,6,npad,0);
% 
% % Compute tau_p from p-pick on
% tauP_h = zeros(1,n);
% Vi     = zeros(1,n);
% Ai     = zeros(1,n);
% 
% %for i = ppxIdx+1:n
% for i = 2:n
%     Vi(i)     = alpha*Vi(i-1) + vel_h(i)^2;
%     Ai(i)     = alpha*Ai(i-1) + acc_h(i)^2;
%     tauP_h(i) = 2*pi*sqrt(Vi(i)/Ai(i));
% end
% 
% % Find max-value of tauP_l in first 0.5 - 4 seconds since p-pick
% startIdx_h = ppxIdx + nsSkip;
% endIdx_h   = ppxIdx + 4*sr;
% 
% [tauP_max_h,idx_h] = max(tauP_h(startIdx_h:endIdx_h));
% m_h                = 7.0*log10(tauP_max_h)+5.9;




%% FIGURE -----------------------------------------------------------------
if (o_plot); 
    hf      = figure(311); clf; whitebg('w')
    hsub(1) = subplot(3,1,1);
    hsub(2) = subplot(3,1,2);
    hsub(3) = subplot(3,1,3);
    
    axes(hsub(1));  hold on
    mm         = 7.14*log10(tau_p) + 5.93;    % Coeffs from Olson and Allen 2005
    ymax       = 1.1*tauPmax;
    l1         = line([tppx tppx], [0 ymax],'Color','r','lineWidth',2*lnWidth);
    [ax,h1,h2] = plotyy(t,tau_p,t,mm);
                 scatter(ax(1),t(idxMax),tauPmax   ,'lineWidth',2,'markerEdgeColor','k')
                 hold(ax(2));  
                 scatter(ax(2),t(idxMax),mm(idxMax),'lineWidth',2,'markerEdgeColor','k')
                 line('parent',ax(2),'xdata',[min(t) max(t)],'ydata',[m m],'color',[0 0.4 0])
                 line('parent',ax(1),'xdata',[tppx+startTime tppx+startTime],'ydata',[0 ymax],'color',[0 0.4 0])
                 line('parent',ax(1),'xdata',[tppx+tauP0     tppx+tauP0]    ,'ydata',[0 ymax],'color',[0 0.4 0])
    set(h2,'lineStyle','-','lineWidth',2*lnWidth,'color',[0 0.4 0])
    set(h1,'lineWidth',2*lnWidth,'color','b')
    set(get(ax(1),'Ylabel'),'String','Tau_p  [s]','fontSize',ftSize)
    set(get(ax(2),'Ylabel'),'String','Magnitude','fontSize',ftSize)
    set(ax(1),'ylim', [0 ymax],'ytick',0:.5:4,'fontSize',ftSize)
    set(ax(2),'ylim', [3 9],'ytick',3:1:9,'fontSize',ftSize)
    set(ax(1),'fontSize',ftSize)
    
    % Tau_p       -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
%     subplot(3,1,1)
%     hold on; grid on;
%     ymax = 1.1*max(abs(tau_p));
%     p1 = plot(t,tau_p, 'k','lineWidth',lnWidth2);
%     l1 = line([tppx tppx], [0 4],'Color','r','lineWidth',2*lnWidth);
%     l2 = line([tppx+startTime tppx+startTime], [0 4],'Color','b','lineWidth',lnWidth);
%     l3 = line([tppx+tauP0   tppx+tauP0],   [0 4],'Color','k','lineWidth',lnWidth);
%     plot(t(idx+ppxIdx+startTime*sr),tauPmax,'or','markerSize',10,'lineWidth',2)
%     ylabel('Tau_p [sec]','fontSize',ftSize)
    
    tStrng = sprintf('Tau_p using alpha= %4.2f and tauP0= %3.1f sec for a M %3.1f record at %3.0f km',alpha,tauP0,m,r);
    title(tStrng,'fontSize',ftSize)

%    set(gca,'ylim',[0 4],'xlim',[min(t) max(t)],'fontSize',ftSize);
    %lgd = legend([p1, p2],'tau_p\_max\_l','tau_p\_max\_h');
    %set(lgd,'Location','NorthWest')
    if nargin>5; set(ax,'xlim',[tppx-ts tppx+te]); end

    % Waveform    -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    axes(hsub(2))
    plot(t,acc,'k','lineWidth',lnWidth)
    hold on; grid on;
    ymax = 1.1*max(abs(acc));
    line([tppx tppx], [-ymax ymax],'Color','r','lineWidth',2*lnWidth)
    ylabel('Acceleration [m/s^2]','fontSize',ftSize)
    set(gca,'fontSize',ftSize)
    set(gca,'ylim',[-ymax ymax])
    set(gca,'xlim',[min(t) max(t)])
    if nargin>5;  set(gca,'xlim',[tppx-ts tppx+te]); end
    
    axes(hsub(3))
    plot(t,vel,'k','lineWidth',lnWidth)
    hold on; grid on;
    ymax = 1.1*max(abs(vel));
    line([tppx tppx], [-ymax ymax],'Color','r','lineWidth',2*lnWidth)
    ylabel('Velocity [m/s]','fontSize',ftSize)
    xlabel('Time [sec]','fontSize',ftSize)
    set(gca,'fontSize',ftSize)
    set(gca,'ylim',[-ymax ymax])
    set(gca,'xlim',[min(t) max(t)])
    if nargin>5;  set(gca,'xlim',[tppx-ts tppx+te]); end

    out.hf = hf;
end