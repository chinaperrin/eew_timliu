function [amax,amaxIdx,cav,noise,sout] = fbank_1trace(s,ppxIdx,sr,tnoise,tmax,ns,fc,fMode,fOrder)
% Returns maximuim amplitudes in all frequency-band specified by fc between
% the pick-time <tpx> and <tpx + snippetLength>, i.e. in increasingly long 
% time windows.

% INPUTS 
% s          = signal/waveform
% sr         = sampling rate
% tnoise     = length of noise window in seconds
% ppxIdx     = index of pick-time
% ns         = number of samples per snippet
% nsnip       = no. of snippets to be evaluated; if empty, uses max. no. of 
%              snippets; if full waveform should be processed as 1 snippet, 
%              choose ns = numel(s(ppxIdx:end))-1 and nsnip=1 as input
%              arguments
% fc         = passband corner frequencies


fprintf(1,'\nNote: This function still uses the principle of \nmeasure_peakAmps_in_expanding_tWindows.m \nwhich makes rounding errors for sr<100sps.\n Change to using the principle of measure_tdpa.m instead. Do it. Now.\n')
nband     = size(fc,1);
nnw       = tnoise*sr;     % No. of samples in noise window
pctl      = [84.1, 97.7];   % Percentiles for measuring noise
o.verbose = false;

% If nsnip is empty, fit as many snippets as possible between p-px and signal end
%if (isempty(nsnip))
%    nsnip = floor((numel(s) - ppxIdx - 1)/ns)-1;
%end


% Check if pre-pick signal is long enough to include <tnoise>
n      = numel(s);
nsnip1 = floor((ppxIdx-1)/ns);         % How many snippets fit between signal start and ppxIdx?
nsnip2 = tnoise*sr/ns;                 % How many snippets fit witin the first 2sec?

if nsnip2<nsnip1; nsnipNoise = nsnip2; % If much more snippets than fit into tnoise are available, use only tnoise
else              nsnipNoise = nsnip1; % If less is available than what would fit into tnoise, use as many as there are
end

% Compute indices of snippet edges
sIdx          = ppxIdx - nsnipNoise*ns;         % Starting index of first snippet
snipEdges     = sIdx-1:ns:n;                    % Edges of all snippets
idxMax        = ppxIdx+tmax*sr;                 % To save time, stop measuring at <tmax> seconds after p-pick
snipEdges     = snipEdges(snipEdges<=idxMax);
snipEdges_cav = snipEdges(2:end) - ppxIdx + 1;      % Snippet edges for cav


% DONT KNOW IF THE +1 IS CORRECT    CHECK   CHECK   CHECK   CHECK   CHECK   

% Initiate arrays
nsnip   = numel(snipEdges)-1;
sout    = cell(nband,1);
amax    = zeros(nband,nsnip);
amaxIdx = zeros(nband,nsnip);
cav     = zeros(nband,nsnip);
noise   = zeros(nband,numel(pctl));


for iband = 1:nband
    
    % * * * * * * * * * * * * * * * * * * * * * *~ * * * * * * * * * * *
    fLo         = fc(iband,1);     % Lower corner frequency (for hipass)
    fUp         = fc(iband,2);     % Upper corner frequency (for lopass)
    sout{iband} = bworth(s,sr,[fLo,fUp],'band',fOrder,fMode);
    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    % Compute maxima in each of the <nInc> time snippets
    for isnip = 1:nsnip
        interval                        = sIdx:snipEdges(isnip+1); % Indices of snippet elements
        S                               = sout{iband}(interval);   % Snippet
        [amax(iband,isnip),amaxIdx_tmp] = max(abs(S));             % Peak amps
        amaxIdx(iband,isnip)            = amaxIdx_tmp + ppxIdx-1;  % Indices of peak amps wrt/ wform start
    end
    
    % Compute CAV and save resp values at each snippet end
    % Note: after snippet-loop, S has it's maximal lenght
    S_cav        = cumsum(abs(S))/sr;
    cav(iband,:) = S_cav(snipEdges_cav);
    
    % Measure noise
    if (ppxIdx>nnw); noiseAmps = sout{iband}(ppxIdx-nnw :ppxIdx-1);
    else             noiseAmps = sout{iband}(1          :ppxIdx-1);
    end
    noise(iband,:) = prctile(noiseAmps,pctl);
end