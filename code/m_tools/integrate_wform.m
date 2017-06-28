function [sout] = integrate_wform(s,sr,ppxIdx,o)
% Integrates time series using cumtrapz, either starting at first sample,
% or at the <ppxIdx>^th sample. In the latter case, all previous values to 
% zero.

intMode = o.intMode;
if ~isfield(o,'ns'); ns = [];
else                 ns = o.ns;
end 

if (nargin<3); 
    intMode = 'allWform';
    fprintf(1,'Neither ppx nor integration mode specified. Integrating all of waveform.\n')
end

dt = 1/sr;
 

% A. integrate only from pick on
% ------------------------------
if strcmp(intMode,'afterPx')
   
    % Option 1: make sure first sample is at zero
    sint               = cumtrapz(s(ppxIdx-1:end))*dt;
    offset             = sint(1);
    sout               = zeros(size(s));
    sout(ppxIdx-1:end) = sint-offset;
    
    % Option 2: remove average amplitude from before pick
    % Option 3: ...
    
    
% B. integrate all waveform
% -------------------------
elseif strcmp(intMode,'allWform')
    sout = cumtrapz(s)*dt;
    %sout = cumsum(s).*dt;


% C. integrate from n samples before ppxIdx
% -----------------------------------------
elseif strcmp(intMode,'nsBeforePx')

    if isempty(ns); fprintf(1,'Need to specify parameter <ns> for this integration mode\n'); end
    % Option 1: make sure first sample is at zero
    sIdx = ppxIdx-ns;
    if sIdx<1; sIdx=1; 
               fprintf(1,'8ung: Wform too short for nsBeforePx. Starting integration at first sample.\n')
    end
    sint = cumtrapz(s(sIdx:end))*dt;
    sout = zeros(size(s));
    %sout(ppxIdx-ns:end) = sint-sint(1);    % Set amplitude at first sample to zero
    sout(ppxIdx-ns:end) = sint-sint(ns);     % Set amplitude at tppx to zero
    
else
    error('Specify integration mode\n')
end