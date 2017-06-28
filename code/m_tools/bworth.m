function [sOut] = bworth(sIn,sr,fCorner,fType,fOrder,fMode)

% sIn       input signal
% sr        sampling rate
% fCorner   upper corner freq fUp for lowpass filter
%           lower corner freq fLo for highpass filter
%           [fLo, fUp]            for bandpass filter
% fMode     'causal' or 'acausal'
%
% menandrin@gmail.com, 140910

fNyq = sr/2;
wNyq = fNyq*2*pi;  % Circular Nyquist frequency

% Filter coefficients
if ( strcmp(fType,'high')  || strcmp(fType,'band') )
    
    wLow  = fCorner(1)*2*pi;                    % Circular corner frequencies
    [B,A] = butter(fOrder,wLow/wNyq,'high');    % Compute coefficients
    
elseif strcmp(fType,'low')
    
    wUp   = fCorner(1)*2*pi;                    % Circular corner frequencies
    [B,A] = butter(fOrder,wUp/wNyq);            % Compute coefficients
end

% Run filter
if     strcmp(fMode,'causal');  sOut = filter  (B,A,sIn);
elseif strcmp(fMode,'acausal'); sOut = filtfilt(B,A,sIn);
end
    

if strcmp(fType,'band')
    
    % Run filter again, this time lowpass, using output from first filter passing as input
    sIn = sOut;

    wUp   = fCorner(2)*2*pi;                    % Circular corner frequencies
    [B,A] = butter(fOrder,wUp/wNyq);            % Compute coefficients

    if     strcmp(fMode,'causal');  sOut = filter  (B,A,sIn);
    elseif strcmp(fMode,'acausal'); sOut = filtfilt(B,A,sIn);
    end
end





%% APPENDIX
% Computation times for various ways to compute the sample-by-sample filter
% operation. Elapsed times correspond to single filter passes (highpass).
% Bottomline: use filter.m and filtfilt.m, they are way faster than anything else.

% % Zero padding
% sIn  = [zeros(fOrder,1); sIn]; 
% ns   = length(sIn);
% sOut = zeros(ns,1);

% Direct form filter
% ixb0 = fOrder+1:-1:1;
% ixa0 = fOrder:-1:1;

% for is = fOrder+1:ns
%     
%     % i. Using flipud                                                    Elapsed time is 0.574510 seconds.
%     % sOut(is) = B*flipud(sIn(is-fOrder:is)) - A(2:end)*flipud(sOut(is-fOrder:is-1));
%     
%     %         % ii. Using precomputed indices                            Elapsed time is 0.385536 seconds.
%     %         ixb      = [is  :-1:is-fOrder];
%     %         ixa      = [is-1:-1:is-fOrder];
%     %         sOut(is) = B*sIn(ixb) - A(2:end)*sOut(ixa);
%     %
%     % iii. Using precomputed indices II                                  Elapsed time is 0.397691 seconds.
% %     ixb2      = is-fOrder-1+ixb0;
% %     ixa2      = is-fOrder-1+ixa0;
% %     sOut(is) = B*sIn(ixb2) - A(2:end)*sOut(ixa2);
% 
%     % iv. Explicit formulation of sum                                    Elapsed time is 0.005247 seconds.
%     %sOut(is) = B(1)*sIn(is) + B(2)*sIn(is-1)  + B(3)*sIn(is-2)  + B(4)*sIn(is-3)  + B(5)*sIn(is-4) ...
%     %                        - A(2)*sOut(is-1) - A(3)*sOut(is-2) - A(4)*sOut(is-3) - A(5)*sOut(is-4);
% 
%     % v. Explicit formulation of sum with loop over fOrder               Elapsed time is 0.028545 seconds.
%     tmp  = 0;
%     A(1) = 0;
%     for io = 1:fOrder+1
%         tmp = tmp + B(io)*sIn(is-io+1) - A(io)*sOut(is-io+1);
%     end
%     sOut(is) = tmp;
%     
%     % vi. Convolution?
% 	% Doesn't work because of recursiveness
%     
%     % Other ways?
% 
% end
% sOut = sOut(fOrder+1:end);

% vii. filter.m                                                            Elapsed time is 0.003784 seconds.
% A(1) = 1;
% sOut = filter(B,A,sIn);


% Plot frequency response
% [y ff]    = freqz(B,A,4096,sr);
% [h6,y,ff] = plot_freqResp(A,B,sr,fUp,figNum,0);