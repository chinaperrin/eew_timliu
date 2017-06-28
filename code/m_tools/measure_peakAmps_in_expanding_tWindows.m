function [amax,nmax,amaxIdx] = measure_peakAmps_in_expanding_tWindows(s,ppxIdx,sr,tnoise,tmax,ns)

% Returns peak absolute amplitudes in increasingly long time windows
% starting at <tnoise> seconds before the p-onset. peak amp vector is then
% split into a pre- and a post-pick part.

% INPUTS 
% s        = signal/waveform
% ppxIdx   = index of start of time windows (e.g. a pick-index)
% sr       = sampling rate
% tnoise   = time interval before pick in which peak amps are measured
% ns       = no. of samples per snippet (--> defines snippet length)

% sIdx                        ppxIdx
%                               |
%                               v
%
%  |        ||        ||       ||        ||        ||        ||        |        
%  ---------  --------  -------  --------  --------  --------  --------  ...
%      n3         n2       n1       s1         s2        s3       s...
%
%  ppxIdx = first sample of first signal window


% Check if pre-pick signal is long enough to include <tnoise>
n      = numel(s);
nsnip1 = floor((ppxIdx-1)/ns);         % How many snippets are available between signal start and ppxIdx?
nsnip2 = floor(tnoise*sr/ns);          % How many snippets fit witin the first <tnoise> sec?

if nsnip2<nsnip1; nsnipNoise = nsnip2; % If much more snippets than fit into tnoise are available, use only tnoise
else              nsnipNoise = nsnip1; % If less is available than what would fit into tnoise, use as many as there are
end

% Compute indices of snippet edges
sIdx      = ppxIdx - nsnipNoise*ns;
snipEdges = sIdx-1:ns:n;
idxMax    = ppxIdx+tmax*sr;                 % To save time, stop measuring at <tmax> seconds after p-pick
snipEdges = snipEdges(snipEdges<=idxMax);

% Initiate arrays
nsnip   = numel(snipEdges)-1;
amax    = zeros(nsnip,1);
amaxIdx = zeros(nsnip,1);

% Compute maxima in each of the <nInc> time snippets
for isnip = 1:nsnip
    interval                     = sIdx:snipEdges(isnip+1); % Indices of snippet elements
    stmp                         = s(interval);             % Signal snippet
    [amax(isnip),amaxIdx(isnip)] = max(abs(stmp));          % Maximum amplitudes
end

% Maxima have been counted wrt/ sIdx, not wrt/ waveform start. In order
% to get amaxIdx rel. to waveform start (index = 1), add sIdx to amaxIdx
amaxIdx = amaxIdx + sIdx-1;

% Split vectors into pre- and post-pick part
nmax    = amax(1:nsnipNoise);
amax    = amax(nsnipNoise+1:end);
amaxIdx = amaxIdx(nsnipNoise+1:end);