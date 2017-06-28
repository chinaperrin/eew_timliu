function [amax,nmax,amaxIdx,sout,tsnippet] = measure_tdpa(s,ppxIdx,sr,tnoise,tmax,snpLength,fbands)
%MEASURE TIME-DEPENDENT (ABSOLUTE) PEAK AMPLITUDES ('tdpa')
%
% Measure tdpa over increasingly long time windows ('snippets'), starting at <tnoise> 
% seconds before the p-onset. The peak amp vector is then split into a pre-
% and a post-pick part. The indices <amaxIdx> are wrt/ to the snippet 
% vector, i.e. if amaxIdx(i)=1 it means that the observed peak amplitude 
% amax(i) occured during the first snippet. 
% If filter corner frequencies are provided in fbands, signal is bandpass
% filtered, using each line of fbands as corner frequencies. If fbands is
% empty, no filters are applied.
%
% INPUTS 
% s             signal time series
% ppxIdx        index of start of time windows (e.g. a pick-index)
% sr            sampling rate in [Hz]
% tnoise,tmax 	peak amps are measured in time interval [-tnoise, tmax] in [sec]
% snpLength     increment by which time-window is increased in [sec]
% fbands        frequency bands in each of which amax is computed
%
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

global fOrder fMode

opt.plotAmps = 0;       % For debugging 


%% Frequency organisation
if   isempty(fbands); nbands          = 1;
                      flg_applyFilter = false;
else                  nbands          = size(fbands,1);
                      flg_applyFilter = true;
end


%% Time organisation
dt = 1/sr;
ns = numel(s); 
ts = (1:ns)'*dt;        
ts = ts-ts(ppxIdx)+dt;  % --> time vector at sampling rate with t=0 at last sample before pick

% truncate signal and time vector to interval [-tnoise, tmax]
isInsideInterval = logical(ts>-tnoise &ts<=tmax);
ts               = ts(isInsideInterval); 
s                = s (isInsideInterval);

% snpLength = 'fullWform' if peak amps should be evaluated over entire signal length
%if ~isnumeric(snpLength); snpLength = ts(end)-dt; end
if strcmp(snpLength,'fullWform'); snpLength = ts(end)-dt; end

tsnippet_noise  = flipud((0:-snpLength:ts(2))'); 
tsnippet_signal = (snpLength:snpLength:ts(end))';
tsnippet        = [tsnippet_noise; tsnippet_signal];
tsnippet        = tsnippet(tsnippet>=-tnoise);
nsnip           = numel(tsnippet);                  % No. of snippets
nsnip_noise     = numel(tsnippet_noise);

if opt.plotAmps
    clf; hold on; grid on
    plot(ts,s,'-xk','lineWidth',2)    
    xlm = get(gca,'xlim');
    xtx = xlm(1):snpLength:xlm(2);
    set(gca,'xtick',xtx)
end



% For each snippet, find signal-vector index <eIdx> of last sample that is
% still inside the isnip^th snippet. Peak amps will then be measured in the 
% interval signal(1:eIdx);
amax    = zeros(nbands,nsnip);
amaxIdx = zeros(nbands,nsnip);
sout    = cell (nbands,1    );

for iband = 1:nbands

    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    if flg_applyFilter; fLo  = fbands(iband,1);    % Lower corner frequency (for hipass)
                        fUp  = fbands(iband,2);    % Upper corner frequency (for lopass)
                        stmp = bworth(s,sr,[fLo,fUp],'band',fOrder,fMode);
    else                stmp = s;
    end
    sout{iband} = stmp;
    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    for isnip = 1:nsnip
        
        eIdx      = find(ts<=tsnippet(isnip),1,'last');
        interval  = 1:eIdx;
        sInterval = sout{iband}(interval);      % Signal snippet
        [val,idx] = max(abs(sInterval));        % Maximum absolute amplitudes
        
        if isempty(val); val=0; idx=0; end
        amax   (iband,isnip) = val;
        amaxIdx(iband,isnip) = idx;
    end
    
    if opt.plotAmps; gcf; hold on;
                     plot(tsnippet, amax,'-dm')
                     plot(tsnippet,-amax,'-dm')
    end    
end


% Separate signal and noise amplitudes
% [amax, amaxIdx, tsnippet, isNoise]
isNoise = tsnippet<=0;
nmax    = amax(:, isNoise);
amax    = amax(:,~isNoise);

% So far, amaxIdx was rel. to start of selected time window, i.e. wrt/ ts(1); 
% shift it to make it rel. to ppxIdx:
amaxIdx = amaxIdx(:,~isNoise)-tnoise*sr;