function [signal,noise,idx] = get_tdpa(trList,mRange,rRange,tRange,band)

% Return peak amplitudes from amax, Pd, Pv & Pa fields of traceList from
% specified magnitude-, distance- and time-ranges. If tRange only includes
% positive times, noise is returned empty. Note that noise is only saved
% for Pd, Pv & Pa, but not for narrow-band peak amplitudes <amax>.
% 
% menandrin@gmail.com, last update: 150509

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
    
    if tRange(1)>=0; returnNoise = false;
    else             returnNoise = true;
    end
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
 
% Make all signal traces have same length:
% (1) if they are longer  than tRange(2), truncate them
% (2) if they are shorter than tRange(2), replicate the last values
nmax = floor(tRange(2)/snpLength);
ntr  = numel(idx);
S    = zeros(ntr,nmax);
 
for itr = 1:ntr
    ns = length(amps{itr});
    if ns>=nmax; S(itr,:) = amps{itr}(1:nmax);                                  % (1)
    else         S(itr,:) = [amps{itr}, repmat(amps{itr}(end),1,nmax-ns)];      % (2)
    end
end

% Write to output structure
signal.amps = S;
signal.x    = (1:length(S(1,:)))';   % x-axis
signal.t    = signal.x*snpLength;




%% Noise 
if returnNoise
    % Read amplitude values
    if     flg_usePd;   amps = trList.dspNoise(idx);
    elseif flg_usePv;   amps = trList.velNoise(idx);
    elseif flg_usePa;   amps = trList.accNoise(idx);
    end
    
    % Make all noise traces have same length:
    % (1) if they are longer  than dt = abs(tRange(1)), truncate them
    % (2) if they are shorter than dt = abs(tRange(1)), replicate the last values
    nmax = floor(abs(tRange(1))/snpLength);
    N    = zeros(ntr,nmax);
    
    for itr = 1:ntr
        nn = length(amps{itr});
        if nn>nmax; N(itr,:) = amps{itr}(1:nmax);
        else        N(itr,:) = [repmat(amps{itr}(1),1,nmax-nn),amps{itr}'];
        end
    end
        
    % Write to output structure
    noise.amps = N;
    noise.x    = [-nmax+1:1:0]';     % x-axis
    noise.t    = noise.x*snpLength;

else
    noise = [];
end