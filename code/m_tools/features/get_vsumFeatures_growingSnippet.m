function [features,iE] = get_vsumFeatures_growingSnippet(acc,vel,dsp,sr,snpLength,nsnp)
%Compute various waveform statistics on increasingly long signal snippets
%of vector sums of the three components X, Y, Z in three unit dimensions
%(acc, vel & dsp)

% nsnp                                  No. of snippets over which statistics are computed
% snpLength                             length of each snippet in sec 
nspersnippet = snpLength*sr;            % No. of samples per snippet
nsnpMax      = floor(numel(acc)/nspersnippet);
if nsnp>nsnpMax; nsnp=nsnpMax; end
iE           = (1:nsnp)'*nspersnippet;  % Indices of End of Waveform snippet


% 1. Get waveform features from single vector ..............................
if ~isempty(acc); [accFeatures] = get_singleTrace_features(acc,iE,sr,nsnp); features.acc = accFeatures; end
if ~isempty(vel); [velFeatures] = get_singleTrace_features(vel,iE,sr,nsnp); features.vel = velFeatures; end
if ~isempty(dsp); [dspFeatures] = get_singleTrace_features(dsp,iE,sr,nsnp); features.dsp = dspFeatures; end


% 2. Get multi-vector waveform features ....................................
% if ~isempty(acc) &~isempty(vel); 
%     taovaS = zeros(nsnp,1);
%     for isnp2 = 1:nsnp
%         taovaS(isnp2) = 2*pi*(max(vel(1:iE(isnp2))) - max(acc(1:iE(isnp2))));
%     end
%     features.taova = taovaS;
% end


function featureList = get_singleTrace_features(s,iE,sr,nsnp)
    
    meanV     = zeros(nsnp,1);
    medianV   = zeros(nsnp,1);
    varV      = zeros(nsnp,1);
    stdV      = zeros(nsnp,1);
    p25V      = zeros(nsnp,1);
    p75V      = zeros(nsnp,1);
    iqrV      = zeros(nsnp,1);
    minimax   = zeros(nsnp,1);
    rmsV      = zeros(nsnp,1);
    kurtosisV = zeros(nsnp,1);
    skewnessV = zeros(nsnp,1);
    zcrV      = zeros(nsnp,1);
    caxV      = zeros(nsnp,1);
    madV      = zeros(nsnp,1);
    f38V      = zeros(nsnp,1);
    k2V       = zeros(nsnp,1);
    qtr       = zeros(nsnp,1);
    
    for isnp = 1:nsnp;
        
        interval = 1:iE(isnp);
        stmp     = s(interval);
        n        = numel(stmp);
        T        = n/sr;
    
        % Mean & median
        meanV  (isnp) = mean  (stmp);
        medianV(isnp) = median(stmp);
        varV(isnp)    = var   (stmp);
        stdV(isnp)    = std   (stmp);  % REDUNDANT
        
        % Percentiles & IQR
        p25V(isnp) = prctile(stmp,25);
        p75V(isnp) = prctile(stmp,75);
        iqrV(isnp) = p75V(isnp)-p25V(isnp);
        
        minimax   (isnp) = max(stmp)-min(stmp);
        rmsV      (isnp) = sqrt(mean(stmp.^2));
        kurtosisV (isnp) = kurtosis(stmp);
        skewnessV (isnp)  = skewness(stmp);
        
        s1          = stmp(2:end);
        s2          = stmp(1:end-1);
        crossingIdx = find(s1.*s2<0);
        zcrV(isnp)  = numel(crossingIdx)/T;
       
        % Cumulative absolute acc/vel/dsp
        %caxV(isnp) = cumsum(abs(stmp))/sr;
        caxV(isnp) = sum(abs(stmp))/sr;
    
        % Mean absolute deviation?
        madV(isnp) = sqrt( 1/(n-1)*sum(stmp-mean(stmp)) );
        
        %38: Z_acc_vector_sum = (np.max(np.abs(acc_vector_sum)) - mean_acc_vector_sum) / var_acc_vector_sum
        f38V(isnp) = (max(abs(stmp)) - mean(stmp))/var(stmp);

        %k2_acc_vector_sum = np.power(skew_acc_vector_sum,2) + np.power(kurtosis_acc_vector_sum,2)
        k2V(isnp)  = skewness(stmp)^2+kurtosis(stmp)^2;
        
        % Quarter-ratio = amps in last quarter vs first quarter of signal
        intstart  = 1:floor(1/4*n); 
        intend    = ceil(3/4*n):n;
        qtr(isnp) = median(abs(s(intend)))/median(abs(s(intstart)));
    end
    
    featureList.mean     = meanV;
    featureList.median   = medianV;
    featureList.var      = varV;
    featureList.std      = stdV;
    featureList.p25      = p25V;
    featureList.p75      = p75V;
    featureList.iqr      = iqrV;
    featureList.minimax  = minimax;
    featureList.rms      = rmsV;
    featureList.kurtosis = kurtosisV;
    featureList.skewness = skewnessV;
    featureList.zcr      = zcrV;
    featureList.cax      = caxV;
    featureList.mad      = madV;
    featureList.f38      = f38V;
    featureList.k2       = k2V;
    featureList.qtr      = qtr;
end
end