function [ytap,t] = make_tapered_plaw_onset(ts,te,dt,a1,a2,tFract,o_plot)
% Create tapered power law signal

%dt = 1e-4;
%ts = 0;
%te = .1;
%a1       = 1;
%a2       = 4;  
%tFract = .5;   % Fraction of duration over which taper goes from 1 to 0
if nargin<7; o.plot=false; end
 
t     = (ts:dt:te)';
ttap  = ts + tFract*(te-ts);
nstap = round(ttap/dt);
%nstap = round((te-ttap)/dt);
y     = a1*t.^a2;                               % Power law
%y1(y1<0) = 0;
ytap  = flipud(taper(flipud(y),1/dt,nstap));    % Tapered power law

if o_plot; figure(44); clf; hold on; grid on; box on;
           plot(t,y)
           plot(t,ytap,'r')
end