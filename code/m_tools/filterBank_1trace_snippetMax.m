function [amax,amaxIdx,noise,sout] = filterBank_1trace_snippetMax(s,sr,nspersnp,nsnp,fc,npad)
% Returns maximuim amplitudes in all frequency-band specified by fc in each snippet

% INPUTS 
% s          = signal/waveform
% sr         = sampling rate
% pxIdx      = index of pick-time
% nspersnp   = number of samples per snippet
% nsnp       = no. of snippets to be evaluated; if empty, uses max. no. of 
%              snippets; if full waveform should be processed as 1 snippet, 
%              choose nspersnp = numel(s(ppxIdx:end))-1 and nsnp=1  
% fc         = passband corner frequencies

realisationType = 8;    % Filter implementation type
o_verbose       = false;
nband           = size(fc,1);

% If nsnp is empty, fit as many snippets as possible between p-px and signal end
if (isempty(nsnp))
    nsnp = floor((numel(s) - pxIdx - 1)/nspersnp);
end

% Find indices of snippet ends
lastIdx = pxIdx+nspersnp:nspersnp:pxIdx+nsnp*nspersnp; % End of snippet wrt/ waveform start

% Initiate array
sout    = cell(nband,1);
amax    = zeros(nband,nsnp);
amaxIdx = zeros(nband,nsnp);
noise   = zeros(nband,1);



for iband = 1:nband
    
    % * * * * * * * * * * * * * * * * * * * * * *~ * * * * * * * * * * *
    fLow        = fc(iband,1);    % Lower corner frequency (for hipass)
    fUp         = fc(iband,2);    % Upper corner frequency (for lopass)
    [sout{iband},~,~] = ...
        butter_pass_tdomain_f(s,fLow,fUp,sr,realisationType,npad,o_verbose);
    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    % Compute maxima in each of the <nInc> time snippets
    
    
    HERE: insert efficient signal split formulation
    
    for isnp = 1:nsnp
        
        % Snippet
        S = sout{iband}(pxIdx:lastIdx(isnp));
        
        % Maximum amplitudes
        [amax(iband,isnp),amaxIdx_tmp] = max(abs(S));
        amaxIdx(iband,isnp)            = amaxIdx_tmp + pxIdx-1; 
    end
    
    % Measure noise
    if (pxIdx>3*sr)
        noise(iband) = max(abs(sout{iband}(pxIdx-3*sr:pxIdx-1))) - abs(mean(sout{iband}(pxIdx-3*sr:pxIdx-1)));
    else
        noise(iband) = max(abs(sout{iband}(1:pxIdx-1)))          - abs(mean(sout{iband}(1:pxIdx-1)));
    end
    
end
