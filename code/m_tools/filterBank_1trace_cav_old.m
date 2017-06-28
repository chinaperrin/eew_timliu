function [amax,amaxIdx,cav,noise,sout] = filterBank_1trace_cav(s,sr,ppxIdx,nspersnp,nsnp,fc,npad,nWindow)
% Returns maximuim amplitudes in all frequency-band specified by fc between
% the pick-time <tpx> and <tpx + snippetLength>, i.e. in increasingly long 
% time windows.

% INPUTS 
% s          = signal/waveform
% sr         = sampling rate
% ppxIdx      = index of pick-time
% nspersnp   = number of samples per snippet
% nsnp       = no. of snippets to be evaluated; if empty, uses max. no. of 
%              snippets; if full waveform should be processed as 1 snippet, 
%              choose nspersnp = numel(s(ppxIdx:end))-1 and nsnp=1  
% fc         = passband corner frequencies
% nWindow    = length of noise window in seconds

realisationType = 8;            % Filter implementation type
o.verbose       = false;
nband           = size(fc,1);
nnw             = nWindow*sr;   % Nr. of samples in noise window

% If nsnp is empty, fit as many snippets as possible between p-px and signal end
if (isempty(nsnp))
    nsnp = floor((numel(s) - ppxIdx - 1)/nspersnp)-1;
end

% Find indices of snippet ends
lastIdx     = ppxIdx+nspersnp:nspersnp:ppxIdx+nsnp*nspersnp; % End of snippet wrt/ waveform start
lastIdx_cav = lastIdx-ppxIdx;                         % End of snippet wrt/ ppxIdx

% Initiate array
sout    = cell(nband,1);
amax    = zeros(nband,nsnp);
amaxIdx = zeros(nband,nsnp);
cav     = zeros(nband,nsnp);
noise   = zeros(nband,1);


for iband = 1:nband
    
    % * * * * * * * * * * * * * * * * * * * * * *~ * * * * * * * * * * *
    fLow        = fc(iband,1);    % Lower corner frequency (for hipass)
    fUp         = fc(iband,2);    % Upper corner frequency (for lopass)
    [sout{iband},~,~] = ...
        butter_pass_tdomain_f(s,fLow,fUp,sr,realisationType,npad,o.verbose);
    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    % Compute maxima in each of the <nInc> time snippets
    for isnp = 1:nsnp
        
        % Snippet
        S = sout{iband}(ppxIdx:lastIdx(isnp));
        
        % Maximum amplitudes
        [amax(iband,isnp),amaxIdx_tmp] = max(abs(S));
        amaxIdx(iband,isnp)            = amaxIdx_tmp + ppxIdx-1; 
    end
    
    % Compute CAV and save resp values at each snippet end
    % Note: after snippet-loop, S has it's maximal lenght
    S_cav        = cumsum(abs(S))/sr;
    cav(iband,:) = S_cav(lastIdx_cav);
    
    % Measure noise
    if (ppxIdx>nnw)
        noise(iband) = max(abs(sout{iband}(ppxIdx-nnw :ppxIdx-1))) - abs(mean(sout{iband}(ppxIdx-nnw :ppxIdx-1)));
    else
        noise(iband) = max(abs(sout{iband}(1:ppxIdx-1)))           - abs(mean(sout{iband}(1:ppxIdx-1)));
    end
end