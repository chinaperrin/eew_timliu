function [W] = get_weightVectors(wNames,zList,k,X)

% k = parameter for kNN
% X = gradient ratio from dlog(Y)/dM = X*dlog(Y)/dlog(R)

global o out minimax rEdges mEdges

% Which weights are to be computed or loaded?
wantsW1 = strcmp('W1',wNames);
wantsW2 = strcmp('W2',wNames);
wantsW3 = strcmp('W3',wNames);

mmin = minimax.mmin;
mmax = minimax.mmax;
rmin = minimax.rmin;
rmax = minimax.rmax;

thousands = linspace(1e3,1e6,1e3);
ftSize    = 15;

r      = zList.hypDist;
rprime = log10(zList.hypDist);
m      = zList.m;
n      = numel(m);

% Weight file name
xstrng             = num2str(X,'%4.2f');
xstrng             = strrep(xstrng,'.','p');
weightFileFullName = strcat([out.outDirFullName,'weight/weight_',num2str(n),'_k',num2str(k),'_X',xstrng,'.mat']);
weightFilePathName = strcat([out.outDirFullName,'weight/']);
hasNewWeights      = false;

if logical(exist(weightFileFullName,'file'))
    load(weightFileFullName)
    hasW1 = isfield(W,'w1');
    hasW2 = isfield(W,'w2');
    hasW3 = isfield(W,'w3');
else
    hasW1 = false;
    hasW2 = false;
    hasW3 = false;
end

% Define data-space cells for W1 & W2
dr     = 1;
rEdges = rmin:dr:rmax;
dm     = .1;
mEdges = mmin:dm:mmax;

W.w0 = ones(n,1);

% ========================================================
% W1: Weight based on data counts in non-overlapping cells
% ========================================================
if (wantsW1 && ~hasW1)
    
    fprintf(1,'\nOutdated. You should use log(R) rather than R.\n')
    pause
    
    fprintf(1,'\nComputing weights based on data counts in non-overlapping cells ...\n')
    fprintf(1,['using dr = ',num2str(dr),' and dm = ',num2str(dm),'...\n'])
    
    hasNewWeights = true;
    
    % W1. Weight based on data point counts in non-overlapping cells
    [W1]     = get_weight_matrix_counts(zList.hypDist,zList.m,o.plotWeightMat);
    %W1full   = ceil(weightFactor*W1);
    %W1naught = ones(size(W1));
    
    % Weight per data point
    wt = zeros(n,1);
    for i = 1:n
        wt(i) = get_weight(r(i),m(i),W1);
    end
    
    % Normalise and save to W
    W.w1 = wt./sum(wt);
end



% ===============================================================
% W2: Weight based on Euclidean distance of all other data points
% ===============================================================
if (wantsW2 && ~hasW2)
    
    fprintf(1,'\nOutdated. You should use log(R) rather than R.\n')
    pause
    
    fudge    = 1e-4;
    fprintf(1,['\nComputing weights based on Euclidean distance of all other data points using fudge=',num2str(fudge),'...\n'])
    fprintf(1,['using dr = ',num2str(dr),' and dm = ',num2str(dm),'...\n'])
    
    hasNewWeights = true;
    
    % W2. Weight based on the Euclidean distance of all other data points
    [W2]     = get_weight_matrix_eucDist(r,m,fudge,o.plotWeightMat);
    %W2full   = ceil(weightFactor*W2);
    %W2naught = ones(size(W2));
    
    % Weight per data point
    wt = zeros(n,1);
    for i = 1:n
        wt(i) = get_weight(r(i),m(i),W2);
    end
    
    % Normalise and save to W
    W.w2 = wt./sum(wt);
end


% ========================================
% W3: Weight based on k nearest neighbours
% ========================================
if (wantsW3 && ~hasW3)
    
    %k = 50;
    fprintf(1,['\nComputing weights based on k nearest neighbours using k=',num2str(k),' ...\n'])
    
    hasNewWeights = true;
    
    % W3. Weight based on k nearest neighbours  .  .  .  .  .  .  .  .  .  .  .
    Dk = zeros(n,1);
    for i = 1:n
        if (ismember(i,[1,thousands])); fprintf(1,['k=',num2str(k),'; X=',num2str(X),'\t',num2str(i),' / ',num2str(n),'\n']); end
        
        dr    = X*(rprime(i) - rprime');
        dm    = m(i) - m';
        D     = sqrt(dr.^2+dm.^2);
        Ds    = sort(D);
        Dk(i) = Ds(k);
    end
    
    % Normalise and save to W
    W.w3  = Dk./sum(Dk);
    W.w3k = k;
    W.w3replicate = round(W.w3/min(W.w3));
    
    if o.debugWeight
        fprintf(1,'Run regression_sensitivity_to_weighting.m')
        pause
    end
    
end

 
if ( o.saveWeights && hasNewWeights)
    fprintf(1,['Saving (overwriting) weight-file ',weightFileFullName,'\n'])
    if ~exist(weightFilePathName,'dir')
        unix(['mkdir -p ',weightFilePathName]);
    end
    save(weightFileFullName,'W')
end

if ~hasNewWeights
    fprintf(1,['\tUsing weights as loaded from ',weightFileFullName,'\n'])
end

if o.plotWeightMat
    h=figure(291); clf;
    whitebg('k')
    subplot(2,1,1); hold on
    scatter(r,W.w3replicate,20,m,'filled')
    xlabel('Hypocentral distance [km]','fontSize',ftSize)
    ylabel('Integer weights','fontSize',ftSize)
    tString = strrep(weightFileFullName,'_','-');
    title(tString,'fontSize',ftSize)
    set(gca,'fontSize',ftSize)
    colorbar
    
    subplot(2,1,2)
    scatter(m,W.w3replicate,20,r,'filled')
    xlabel('Magnitude','fontSize',ftSize)
    ylabel('Integer weights','fontSize',ftSize)
    set(gca,'fontSize',ftSize)
    colorbar 
    
    if (o.printWeightMat)
        fprintf(1,'Printing weight vector and figure to file...\n')
        figFileName = strrep(weightFileFullName,'.mat','');
        print('-dpng',figFileName)
        saveas(h,figFileName,'fig') 
        fprintf(1,' done.\n')
    end
end



%% ======
% OTHERS
% ======

% Yavors suggestion:
% estimate log10(R)- and M-density independently with ksmooth
% for each data point, compute pr & pm
% ptot = pr*pm;
% w = 1/ptot;

% % Weighting using Gaussian mixture model
% R = zList.hypDist;
% M = zList.m;
% %options = statset('Display','final','maxIter',1e4);
% options = statset('Display','iter','maxIter',1e4);
% diary('log_llh.txt')
% obj = gmdistribution.fit([log10(R), M],40,'Options',options);
% diary off
% p   = pdf(obj,[log10(R),M]);
% pct = prctile(p,0.05);
% p(p<pct) = pct;
% w = 1./p;
% wn = w./sum(w);
%
%
% % 2D Kernel density estimation
% [bandwidth,p,X,Y] = kde2d([log10(R),M],2^7,[log10(0.1) log10(rmax)],[mmin mmax]);
% figure(222); clf;
% surf(X,Y,p,'LineStyle','none')
% colormap copper
% hold on, alpha(.8)
% plot(log10(R),M,'w.','MarkerSize',5)
%
% pct = prctile(p,0.05);
% p(p<pct) = pct;
% w = 1./p;
% wn = w./sum(w);
%
%
% figure(299); clf; hold on;
% scatter(R,M,15,p,'filled')
% colorbar;
%
% hold off
% scatter(R,M,15,wn,'filled')
% scatter(R,M,15,w,'filled')
% scatter(R,M,15,log10(w),'filled')
% colorbar
