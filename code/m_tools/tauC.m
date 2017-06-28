function [out] = tauC(traceFullName,trList,tau0,o_plot,ts,te)

% Computes tau_c as defined in Wu and Kanamori 2005; Applies highpass 
% filter with fc = 0.075Hz. Params for m-estimations can be found in 
% Wu et al., 2007, GJI: If specified and available, <ts> and <te> give the 
% window that is plotted.
%
% menandrin@gmail.com, 130807

% Notes:
% - Wu and Kanamori 2005 do not lowpass filter waveforms, they only apply a
%   highpass filter with fc = 0.075Hz to stabilise the double integration
% -------------------------------------------------------------------------

global vs


%% 1. Preparations
lnWidth   = 1;
ftSize    = 15;

fOrder  = 4;
fMode   = 'causal';
fCorner = 0.075;

if (nargin<4); o_plot = false; end

% Find entry in traceList and read meta-data
idx       = find(strcmp(traceFullName,trList.fullName));
m         = trList.m(idx);
r         = trList.hypDist(idx);
orntCode  = trList.orntCode{idx};

% Avoid using the S-phase: 
if r<vs*tau0; tau0 = r/vs;
              fprintf(1,'tauC window truncated to avoid using S-phase\n')
end

if (~strcmp(orntCode,'Z'))
    fprintf(1,'Warning: not a vertical component! This is not how tau_c is defined ...\n')
end

% Read and pre-process waveforms (Acc/Vel/Dsp)
[S,meta] = read_any_trace(trList.fullName{idx},trList,1);
% vel      = S.vel;
% dsp      = S.dsp;
vel      = bworth(S.vel,meta.sr,fCorner,'high',fOrder,fMode);   % Highpass filtered after integration
dsp      = bworth(S.dsp,meta.sr,fCorner,'high',fOrder,fMode);
t        = meta.t;
ppxIdx   = meta.ppxIdx;
tppx     = t(ppxIdx);
n        = length(vel);
sr       = meta.sr;
dt       = 1/sr;
ns       = round(tau0*sr)-1;     % No. of samples for tauC estimation

% Compute tauC
interval = [ppxIdx ppxIdx+ns];
r_1      = cumtrapz(vel(interval).^2)*dt;
r_2      = cumtrapz(dsp(interval).^2)*dt;
r        = r_1(end)/r_2(end);
tau_c    = 2*pi/sqrt(r);

% Compute magnitude with coefficients from Wu and Kanamori, 2005, BSSA (only 12 events) 
m_wu05 = 4.525*log10(tau_c) + 5.036;

% Arrange output
out.m_wu05 = m_wu05;
out.tau_c  = tau_c;




%% Plot tauC    -----------------------------------------------------------
if (o_plot);
    
    % Compute tau_c as a time series (rather than a scalar)
    tau_c2       = zeros(n,1);
    m_est        = zeros(n,1);
    isIncomplete = false(n,1);
    
    for i = 2:n                    % Start computing tauC at first samples, even if not <ns> samples are available;
        sIdx      = i - ns + 1;    % just use all available samples
        eIdx      = i;
        
        if sIdx<1 
            sIdx            = 1;
            isIncomplete(i) = true;
        end
        
        r_1       = cumtrapz(vel(sIdx:eIdx).^2)*dt;
        r_2       = cumtrapz(dsp(sIdx:eIdx).^2)*dt;
        r         = r_1(end)/r_2(end);
        tau_c2(i) = 2*pi/sqrt(r);
        m_est(i)  = 4.525*log10(tau_c2(i)) + 5.036;     % Wu and Kanamori, 2005, BSSA (only 12 events)
    end

    
    % Plot
    hfig1   = figure(211); clf; whitebg('w');
    hsub(1) = subplot(3,1,1);
    hsub(2) = subplot(3,1,2);
    hsub(3) = subplot(3,1,3);
    
    % Tau_c       -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    axes(hsub(1))
    hold on; grid on;
    
    ymax       = 1.1*max(tau_c2);
    l1         = line([tppx tppx], [0 ymax],'Color','r','lineWidth',2*lnWidth);
    [ax,h1,h2] = plotyy(t,tau_c2,t,m_est);
    plot(t(isIncomplete),tau_c2(isIncomplete),'.k')
    
    set(ax,{'ycolor'},{'b'; [0 0.4 0]})  % Left color red, right color blue...
    set(h1,'color','b')
    set(h2,'color',[0 0.4 0])

    % Start and end indices of plotting window
    if nargin>4;    sIdx = ppxIdx-ts*sr;
                    eIdx = ppxIdx+te*sr;
    else            sIdx = 1;
                    eIdx = n;
    end
    if sIdx<1; sIdx=1; end  % Make sure indices are not outside waveform size
    if eIdx>n; eIdx=n; end
    
    % Plot time and tauC value that was selected with specified tau0
    estIdx = ppxIdx + ns;   
                 scatter(ax(1),t(estIdx),tau_c2(estIdx),'lineWidth',2,'markerEdgeColor','k')
    hold(ax(2)); scatter(ax(2),t(estIdx),m_est (estIdx),'lineWidth',2,'markerEdgeColor','k')
                 line('parent',ax(2),'xdata',[min(t) max(t)],'ydata',[m m],'color',[0 0.4 0])
    
    set(h2,'lineStyle','-','lineWidth',2*lnWidth,'color',[0 0.4 0])
    set(h1,'lineWidth',2*lnWidth,'color','b')
    set(get(ax(1),'Ylabel'),'String','Tau_c  [s]','fontSize',ftSize)
    set(get(ax(2),'Ylabel'),'String','Magnitude','fontSize',ftSize)
    set(ax(2),'ylim', [3 9],'ytick',3:1:9,'fontSize',ftSize)
    set(ax(1),'fontSize',ftSize)
    
    tStrng = sprintf('Tau_c using %3.1fsec for a M%3.1f record at %3.0fkm',tau0,m,r);
    title(tStrng,'fontSize',ftSize)
    set(ax,'xlim',[min(t) max(t)])
    if nargin>4; set(ax,'xlim',[tppx-ts tppx+te]); end
    
    
    
    % Waveform    -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    axes(hsub(2))
    plot(t,vel,'k','lineWidth',lnWidth)
    hold on; grid on;
    ymax = 1.1*max(abs(vel(sIdx:eIdx)));
    line([tppx tppx], [-ymax ymax],'Color','r','lineWidth',2*lnWidth)
    ylabel('Velocity [m/s]','fontSize',ftSize)
    set(gca,'fontSize',ftSize)
    set(gca,'ylim',[-ymax ymax])
    set(gca,'xlim',[min(t) max(t)])
    if nargin>4; set(gca,'xlim',[tppx-ts tppx+te]); end
    
    axes(hsub(3))
    plot(t,dsp,'k','lineWidth',lnWidth)
    hold on; grid on;
    ymax = 1.1*max(abs(dsp(sIdx:eIdx)));
    line([tppx tppx], [-ymax ymax],'Color','r','lineWidth',2*lnWidth)
    ylabel('Displacement [m]','fontSize',ftSize)
    xlabel('Time [sec]','fontSize',ftSize)
    set(gca,'fontSize',ftSize)
    set(gca,'ylim',[-ymax ymax])
    set(gca,'xlim',[min(t) max(t)])
    if nargin>4; set(gca,'xlim',[tppx-ts tppx+te]); end
    
    out.hf = hfig1;
end