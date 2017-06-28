function [amax,nmax,amaxIdx] = measure_tdpa(s,ppxIdx,sr,tnoise,tmax,snpLength)
%
% MEASURE TIME-DEPENDENT (ABSOLUTE) PEAK AMPLITUDES ('tdpa')
%
% Measure tdpa over increasingly long time windows ('snippets'), starting at <tnoise> 
% seconds before the p-onset. The peak amp vector is then split into a pre-
% and a post-pick part. The indices <amaxIdx> are wrt/ to the snippet 
% vector, i.e. if amaxIdx(i)=1 it means that the observed peak amplitude 
% amax(i) occured during the first snippet. 
%
% INPUTS 
% s             signal time series
% ppxIdx        index of start of time windows (e.g. a pick-index)
% sr            sampling rate in [Hz]
% tnoise,tmax 	peak amps are measured in time interval [-tnoise, tmax] in [sec]
% snpLength     increment by which time-window is increased in [sec]
%
% ISSUE
% Not sure if amaxIdx is correct. I think it is relative to the sInterval
% vector, not rel. to the tsnippet vector
%
%                            ppxIdx
%                               |
%                               v
%
%  |        ||        ||       ||        ||        ||        ||        |        
%  ---------  --------  -------  --------  --------  --------  --------  ...
%      n3         n2       n1       s1         s2        s3       s...
%
%  ppxIdx = first sample of first signal window; by definition, ppxIdx is
%           the index of the first sample that is above the noise


opt.plotAmps = 0;       % For debugging

dt = 1/sr;
ns = numel(s); 
ts = (1:ns)'*dt;
ts = ts-ts(ppxIdx)+dt;


% truncate signal to be no longer than [-tnoise, tmax]
isInsideInterval = logical(ts>-tnoise &ts<=tmax);
ts               = ts(isInsideInterval); 
s                = s (isInsideInterval);

% snpLength = 'fullWform' if peak amps should be evaluated over entire
% signal length
if ~isnumeric(snpLength); snpLength = ts(end)-dt; end
%if ~isnumeric(snpLength); snpLength = floor(ts(end)); end

if opt.plotAmps
    clf; hold on; grid on
    plot(ts,s,'-xk','lineWidth',2)    
    xlm = get(gca,'xlim');
    xtx = xlm(1):snpLength:xlm(2);
    set(gca,'xtick',xtx)
end


tsnippet_noise  = flipud((0:-snpLength:ts(1))');
tsnippet_signal = (snpLength:snpLength:ts(end))';
tsnippet        = [tsnippet_noise; tsnippet_signal];
tsnippet        = tsnippet(tsnippet>=-tnoise);
nsnip           = numel(tsnippet);                  % No. of snippets


% For each snippet, find signal-vector index <eIdx> of last sample that is
% still inside the isnip^th snippet. Peak amps will then be measured in the 
% interval signal(1:eIdx);
amax    = zeros(nsnip,1);
amaxIdx = zeros(nsnip,1);

for isnip = 1:nsnip
    
    eIdx      = find(ts<=tsnippet(isnip),1,'last');
    interval  = 1:eIdx;
    sInterval = s(interval);                % Signal snippet
    [val,idx] = max(abs(sInterval));        % Maximum absolute amplitudes
    
    if isempty(val); val=0; idx=0; end
    amax   (isnip) = val;
    amaxIdx(isnip) = idx;
end

if opt.plotAmps
    gcf; hold on;
    plot(tsnippet, amax,'-dm')
    plot(tsnippet,-amax,'-dm')
end

% Separate signal and noise amplitudes
% [amax, amaxIdx, tsnippet, isNoise]
isNoise     = tsnippet<=0;
nsnip_noise = sum(isNoise);
nmax        = amax( isNoise);
amax        = amax(~isNoise);

% So far, amaxIdx was rel. to signal start; shift it to make it rel. to 
% first snippet after p-onset. amaxIdx<=0 means that the observed amplitude
% comes from before the pick. This happens, e.g. if sr=80 and snpLength=.01
% i.e. when the first snippet does not contain a sample from after the
% pick.
amaxIdx     = amaxIdx(~isNoise)-nsnip_noise;