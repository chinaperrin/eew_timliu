function pEx = get_propOfExceedence_ecdf(vals,thresholds)
%PROBABILITY OF EXCEEDENCE USING ECDF
% 
% Copmutes probability of exceedence <pEx> of the levels in <thresholds>, 
% given a vector of realisations of quantity <vals>, using the empirical 
% cumulative distribution function

[ecp,valVect] = ecdf(vals);                                                     % Compute empirical CDF

nthresh = numel(thresholds);
pEx     = zeros(nthresh,1);
for ithresh = 1:nthresh                                                         % For each threshold ...
    
    [~,imin] = min(abs(valVect-thresholds(ithresh)));                           % Find closest match in valVect
    if     imin==1;              pEx(ithresh) = 1;                              % If all values are above threshold, pEX=1
    elseif imin==numel(valVect); pEx(ithresh) = 0;                              % If all values are below threshold, pEX=1
    else
        if valVect(imin)<thresholds(ithresh); imin = [imin  ; imin+1];          % Find the two bracketing values in valVect 
        else                                  imin = [imin-1; imin  ];          % and interpolate to the cumulative prob.
        end                                                                     % at the threshold level
        ecpPrime     = interp1(valVect(imin), ecp(imin),thresholds(ithresh));
        pEx(ithresh) = 1-ecpPrime;
    end
end


% TEST
% For normally distributed samples, 84% of samples lie above y = mu-1*sigma
% mu    = 10;
% sigma = 2;
% nrnd  = 1e6;
% R     = normrnd(mu,sigma,nrnd,1);
% sum(R>(mu-sigma))/nrnd*100;
% 
% thresholds = mu-sigma;
% pEx        = get_propOfExceedence_ecdf(R,thresholds);
% pEx        = pEx*100;