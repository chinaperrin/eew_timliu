function [h1,h2,h3,h4] = reproduce_fbout_1trace(tr_idx,trList,snpLength)

addpath(genpath('../'))

if ~exist('sWindow','var'); sWindow = 1; end
if ~exist('nWindow','var'); nWindow = 1; end

o.verbose = 0;

tFrac = 0.1;
T     = 30;
cmap  = flipud(copper(9));

fOrder = 4;

fc = [24.0     48.0;        % Passbands
      12.0     24.0;
       6.0     12.0;
       3.0      6.0;
       1.5      3.0;
       0.75     1.5;
       0.375    0.75;
       0.1875   0.375;
       0.09375  0.1875];

vp     = 6;              % constant p-phase speed
vs     = 3.4;            % constant s-phase speed
sdelay = (1/vs-1/vp);    % Delay of s-phase wrt/ p-phase per km [s/km]

traceFullName = trList.fullName{tr_idx};
slash_idx     = regexp(traceFullName,'/');
traceName     = traceFullName(slash_idx(end)+1:end);
ornt          = trList.orntCode{tr_idx};
        
        
% Read trace
[S,meta] = read_any_trace(traceFullName,trList,1);
s        = S.vel;
sr       = meta.sr;
ppxIdx   = meta.ppxIdx;
t        = meta.t;

tppx      = t(ppxIdx);
hD        = trList.hypDist(tr_idx);
dtsp      = hD*sdelay;
tspx      = tppx + dtsp;

nspersnp = round(snpLength*sr);

[~,noiseAmps,~,~] = get_snr(s,ppxIdx,sr,sWindow,nWindow);
noise             = prctile(noiseAmps,84.1);

% Pass it through filter bank
[amax,amaxIdx,~,sbNoise,sout] = fbank_1trace(s,ppxIdx,sr,0,1e4,nspersnp,fc,'causal',fOrder);

sIdx        = round(ppxIdx - T*sr*tFrac);      % t-index of window-start
eIdx        = round(ppxIdx + T*sr*(1-tFrac));  % t-index of window-end
titleString = strcat(['Event ',traceName]);

%[h1] = plot_wform(s   ,t,tppx,tspx,591);
h1 = [];
[h2] = plot_wform(S.raw,t,tppx,tspx,592);
if ( strcmp(meta.origin,'scsn') || strcmp(meta.origin,'scsnPx') )
	plot_in_sac(traceFullName,'off')
end

% Plot traces in all passbands
[h3] = plot_nbWforms(s,sout,amaxIdx,t,tppx,tspx,noise,sbNoise,sIdx,eIdx,fc,titleString,tFrac,T);
       

% Plot max-amplitude spectrum
fprintf(1,'plot_nbPGV.m does no longer work. Update it.\n')
%[h4] = plot_nbPGV(amax,fc,sbNoise,[],snpLength,titleString,o.verbose);
h4 = [];

hold on
%plot(amax_obs,'xr','MarkerSize',15) --> for verification