function [signal,noise,idx] = get_time_dependent_peak_amps(trList,mRange,rRange,band,pctl,yscale)

global snpLength


% Which amplitudes to plot?
if isnumeric(band)
    flg_useAmax = true;                 % Use nbPGV as stored in trList.amax
    returnNoise = false;
else
    flg_useAmax = false;
    flg_usePd   = strcmp(band,'Pd');    % Use broad band peak amps
    flg_usePv   = strcmp(band,'Pv');
    flg_usePa   = strcmp(band,'Pa');
    returnNoise = true;
end
    
m   = trList.m;
r   = trList.hypDist;
idx = find( m>=mRange(1) & m<mRange(2) & r>=rRange(1) & r<rRange(2) );


%% Signal
% Read amplitude values
if     flg_useAmax; amps = cellfun(@(x) x(band,:), trList.amax(idx),'uniformOutput',0);
elseif flg_usePd;   amps = trList.pd(idx);
elseif flg_usePv;   amps = trList.pv(idx);
elseif flg_usePa;   amps = trList.pa(idx);
end

% Fill in zeros until all time series are as long as longest one
nmax = max(cellfun(@(x) numel(x), amps));
ntr  = numel(idx);
A    = zeros(ntr,nmax);

for itr = 1:ntr
    na       = length(amps{itr});
    A(itr,:) = [amps{itr}, repmat(amps{itr}(end),1,nmax-na)];
end

if strcmp(yscale,'log'); A = log10(A); end

% Lower, middle and upper percentiles
pl = prctile(A,pctl(1))';
pm = prctile(A,pctl(2))';
pu = prctile(A,pctl(3))';

% Write to output structure
signal.amps = A;
signal.pl   = pl;
signal.pm   = pm;
signal.pu   = pu;
signal.x    = (1:length(pm))';   % x-axis
signal.t    = signal.x*snpLength;




%% Noise
if returnNoise
    % Read amplitude values
    if     flg_usePd;   amps = trList.dspNoise(idx);
    elseif flg_usePv;   amps = trList.velNoise(idx);
    elseif flg_usePa;   amps = trList.accNoise(idx);
    end
    
    % Fill in zeros at start of noise vectors until all time series are as 
    % long as longest one (which has been defined with <tnoise> in the
    % fbank_wformProc scripts.
    nmax = max(cellfun(@(x) numel(x), amps));
    N    = zeros(ntr,nmax);
    
    for itr = 1:ntr
        nn     = length(amps{itr});
        N(itr,:) = [repmat(amps{itr}(1),1,nmax-nn),amps{itr}'];
    end
    
    if strcmp(yscale,'log'); N = log10(N); end
    
    % Lower, middle and upper percentiles
    pl = prctile(N,pctl(1))';
    pm = prctile(N,pctl(2))';
    pu = prctile(N,pctl(3))';
    
    % Write to output structure
    noise.amps = N;
    noise.pl   = pl;
    noise.pm   = pm;
    noise.pu   = pu;
    noise.x   = [-nmax+1:1:0]';     % x-axis
    noise.t   = noise.x*snpLength;

else
    noise = [];
end