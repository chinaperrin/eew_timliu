function [hf,outMat] = cf_peak_amps(trList,mRanges,rRange,tRange,fig)

global fc ftSize iN snpLength colours

%% OVERVIEW
% 1. Preparations
% 2. Reading peak amplitudes from traceList
% 3. Plot PGD(t) curves
%    . PLOT MEDIANS AND RANGES
%    . ADD NOISE PGD(t) CURVES
%    . FORMAT FIGURE
%    . PLOT THEORETICAL S-ARRIVALS
%    . LEGEND
%    . PLOT SLOPE SCALE (outdated)
% 4. KS test
% 5. Plot observed slopes (outdated)
% 6. Print Figure
% 7. Save Output Matrix
% A. Appendix



%% %%%%%%%%%%%%%%%%%
% 1. Preparations %%
%%%%%%%%%%%%%%%%%%%%
fprintf(1,'hardcoded params (should be commented out):\n')
snpLength=.01
ftSize=15


% Quick fix that stops script from crashing when some magnitude bins are empty
nr = size(mRanges,1);
m  = trList.m;
r  = trList.hypDist;

% Exclude mRanges with no data
for ir =1:nr; n(ir) = sum(m>=mRanges(ir,1) &m<mRanges(ir,2)); end
mRanges = mRanges(n>0,:);
nr      = size(mRanges,1);

if strcmp(trList.orntCode{1},'Z'); ornt = 'Vertical';
else                               ornt = 'Horizonal';
end 

if ~isfield(fig,'snrMin'); fig.snrMin=0; end

band   = fig.band;
pctl   = fig.pctl;
plevel = fig.ksPlevel;  % Significance level for estimating separation time tsep, 
                        % = time at which amps in two mBands become significantly
                        %   different, according to KS test.
                
tmin = tRange(1);
tmax = tRange(2);
if tmin==0 & strcmp(fig.xscale,'log'); tmin = .8*snpLength; tRange(1)=tmin; end


% Find out whether it's broadband, strong motion or both
hasH = numel(find(strcmp('H',trList.instrCode)))>0;
hasL = numel(find(strcmp('L',trList.instrCode)))>0;
hasN = numel(find(strcmp('N',trList.instrCode)))>0;
hasG = numel(find(strcmp('G',trList.instrCode)))>0;
hasP = numel(find(strcmp('P',trList.instrCode)))>0;

if (hasH  && ~hasL && ~hasN && ~hasG ); titleAppendix = 'BB'     ; end
if (~hasH && (hasL || hasN || hasG ) ); titleAppendix = 'SM'     ; end
if (hasH  && (hasL || hasN || hasG ) ); titleAppendix = 'SMandBB'; end
if hasP;                                titleAppendix = 'Geophones'; end


% Prepare figure ..........................................................
hf=figure(918); clf; hold on; grid on; box on; whitebg('w')
set(0,'DefaultFigureWindowStyle','docked')

if fig.plotKS;
    % Option 1: plot all KS tests in individual subplot, add subplot for durations
    % Option 2: plot all KS tests in same subplot,       add subplot for durations
    if fig.multipleKSplots; nlines = 4 + nr;
    else                    nlines = 4 + 2;
    end
    hsp = subplot(nlines,6,[1:4,7:10,13:16,19:22]); hold on; grid on; box on;
end

if fig.plotObsSlopes; fprintf(1,'Plotting slope option is outdated. check.\n')
                      hsp(1) = subplot(2,1,1); hold on; grid on; box on;
end





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Reading peak amplitudes from traceList %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'Extracting signal amps from trList .. ')
signal = cell(nr,1);
noise  = cell(nr,1);
idx    = cell(nr,1); 
bt     = cell(nr,1);

% Save station coordinates and magnitudes for plotting map
rec.lat = [];
rec.lon = [];
rec.m   = [];
eqk.lat = [];
eqk.lon = [];
eqk.m   = [];

% Longest pd0 vector; needed as max-dimension for perturbed KS test 
ntmax = max(cellfun(@(x) size(x,2), trList.var4));
pdKS  = cell(nr,1);

for ir = 1:nr
    fprintf(1,sprintf('%i/%i..',ir,nr))
    % Peak amps
    [signal{ir},noise{ir},idx{ir}] = get_tdpa(trList,mRanges(ir,:),rRange,tRange,band);
      
    % Bin titles
    mbin      = floor(mRanges(ir,1));
    remainder = rem(mRanges(ir,1),1);
    if remainder==.5; annotTxt = 'high';
    else              annotTxt = 'low ';
    end
    bt{ir} = sprintf('M%i_{%s}',mbin,annotTxt);
    
    % Prepare KS-test with pick pick-uncertainties
    if fig.considerPxDist
        
        % TraceList for individual magnitude bin
        tmpList  = trList.selectSubList(idx{ir});
        ntmp     = numel(idx{ir});                      % No. of records in bin
        %pdRepMat = [];

        nrnd     = trList.prop.recompPd.params.nrnd;    % No. of picker runs with randmoised starting points
        ntrRep   = ntmp*nrnd;                           % No. of replications
        pdRepMat = zeros(ntrRep,ntmax);                 % Matrix to store Pd-curves from all replications
        ile      = 0;
        ct99999  = 0;
        
        % Construct matrix that has pd-replication for each pick
        for itmp = 1:ntmp
            
            % Pick indices that were found 
            pxSet    = tmpList.var1{itmp}.ppxIdxMat;    % List of picks for this record
            pd0      = tmpList.var4{itmp};              % Pd curve from earliest pick
            uniquePx = unique(pxSet);
            uniquePx = uniquePx(uniquePx~=99999);
            minPx    = min(uniquePx);
            
            ct99999 = ct99999+sum(pxSet==99999);
            
            for ipx = 1:numel(uniquePx)
                ppxIdx = uniquePx(ipx);         % Pick index
                nrep   = sum(pxSet==ppxIdx);    % How often does it occur
                %nrep   = 1;
                ils    = ile+1; 
                ile    = ils+nrep-1;
                dpx    = ppxIdx-minPx;          % Shift rel. to pd0-curve
                is     = 1+dpx; 
                pdShift  = pd0(is:end)-pd0(is);                              % Pd curve starting at ppxIdx
                nshort   = ntmax-numel(pdShift);                     % How many samples shorter than ntmax?
                pdTmp    = [pdShift, repmat(pdShift(end),1,nshort)];  % Make pd time series have length ntmax
                pdRepMat(ils:ile,:) = repmat(pdTmp,nrep,1);
                %ct99999 = ct99999+sum()
            end
            1+1;
        end
        %pdRepMat(:,1:10); clf; proper_hist(pdRepMat)
        %lgc      = (pdRepMat(:,5)~=0);
        %pdKS{ir} = pdRepMat(lgc,:);
        pdKS{ir} = pdRepMat(1:ile,:);
    end
end
ntr = cellfun(@(x) size(x,1), idx);             % No. of traces in each bin
fprintf(1,'done.\n')



%% %%%%%%%%%%%%%%%%%%%%%%%
%  3. Plot PGD(t) curves %
%%%%%%%%%%%%%%%%%%%%%%%%%%

%_ PLOT MEDIANS AND RANGES   . . . . . . . . . . . . . . . . . . . . . . . .
ha     = zeros(nr,1);  %ha      = cell(nr,1);
amps   = cell(nr,1);
tvect  = cell(nr,1);
outMat = cell(nr,1);

% Make sure same m-range has same color
%colours  = {'k','r','c','y','m','b',[0 .4 0],[255/256 165/256 0],[238/256 130/256 238/256],'k','m','c','g','y','b','r'};
colours = fig.colours;
nc      = numel(colours);
%nskipbin = round(2*(6.5-mRanges(1,1)));
%colours  = colours(1+nskipbin:nc);

% Manually plot hypocentral distance ECDF
neverTrue=false;
if neverTrue

    figure(1); clf; hold on; grid on
    for ir = 1:nr; fprintf(1,sprintf('%i/%i\n',ir,nr))
                   [f,x] = ecdf(trList.hypDist(idx{ir}));
                   plot(x,f,'color',colours{ir},'lineWidth',2)
    end
    
    set(gca,'fontSize',ftSize,'YAxisLocation','right','ylim',[0 1],'xlim',rRange)
    ylabel('\Phi (r \leq r_0) ','fontSize',ftSize)
    xlabel('Hypocentral distance [km]','fontSize',ftSize)    
    print('-depsc2','~/Documents/eq/topics/onset/fig/distance_ecdf.eps')
end

for ir = 1:nr
    
    fprintf(1,sprintf('%i/%i\n',ir,nr))
    
    % Unpack signal data
    ts    = signal{ir}.t';
    sAmps = signal{ir}.amps;
    
    nzeros = sum((sAmps(:)==0));
    nnans  = sum(isnan(sAmps(:)));
    if (nzeros~=0 |nnans~=0); fprintf(1,sprintf('8UNG: %i entries with amp=0 and %i entries with amp=Nan\n',nzeros,nnans)); end
    
    % Option 1: wo/ excluding amp-entries which are ==0
    psl = prctile(sAmps,pctl(1));
    psm = prctile(sAmps,pctl(2));
    psu = prctile(sAmps,pctl(3));
    
    % Option 2: Exclude amp-entries which are ==0
    % Moved to APPENDIX
    
    % y-range ("error bars") on lin-scale: in shadedErrorBar.m the error 
    % bars are not to be given in absolute amplitudes, but as differences 
    % wrt/ the middle line (median) --> dps[ul]
    % see m_jumble/understand_shadedErrorBar.m
    dpsu   = psu - psm;        
    dpsl   = psm - psl;      
    sRange = [dpsu; dpsl];
    
    % y-range ("error bars") on log-scale
    pslLog    = log10(psl);
    psmLog    = log10(psm);
    psuLog    = log10(psu);
    dpsuLog   = psuLog - psmLog;
    dpslLog   = psmLog - pslLog;
    sRangeLog = [dpsuLog; dpslLog];
        
    % Unpack noise data
    if ~isempty(noise{1})
        
        tn  = noise{ir}.t';
        pnl = prctile(noise{ir}.amps,pctl(1));
        pnm = prctile(noise{ir}.amps,pctl(2));
        pnu = prctile(noise{ir}.amps,pctl(3));
        
        % y-range ("error bars") on lin-scale
        dpnu   = pnu - pnm;       % Error bars do not have to be given in absolute sense...
        dpnl   = pnm - pnl;       % ... but with respect to middle line
        nRange = [dpnu; dpnl];
        
        % y-range ("error bars") on log-scale
        dpnuLog   = log10(pnu) - log10(pnm);
        dpnlLog   = log10(pnm) - log10(pnl);
        nRangeLog = [dpnuLog; dpnlLog];
        
        % Concatenate noise and signal y-ranges
        t          = [tn        , ts        ];
        pRange     = [nRange    , sRange    ];
        pRangeLog  = [nRangeLog , sRangeLog ];
        pMedian    = [pnm       , psm       ];
        pMedianLog = [log10(pnm), log10(psm)];

        % Save concatenated peak amps in structure
        amps{ir}  = [noise{ir}.amps, signal{ir}.amps];
        tvect{ir} = t;
    else 
        
        % Use signal values only
        t          = ts;
        pRange     = sRange;
        pRangeLog  = sRangeLog; 
        pMedian    = psm;
        pMedianLog = log10(psm);
        
        % Save concatenated peak amps in structure
        amps{ir}  = signal{ir}.amps;
        tvect{ir} = t;
    end 
    
    % Plot percentiles ....................................................
    % 1. Log-log scale  . . . . . . . . . . . . . . . . . . . . . . . . . . 
    if ( strcmp(fig.yscale,'log') & strcmp(fig.xscale,'log') )
        if fig.plotIndividualCurves; plot(t,log10(sAmps)); end
        
                 plot(t,pMedianLog          ,'lineStyle','-','color','w','lineWidth',4);
        ha(ir) = plot(t,pMedianLog          ,'lineStyle','-','color',colours{ir},'lineWidth',2);
        if fig.plotRange
            plot(t,pMedianLog+pRangeLog(1,:),'lineStyle',':','color',colours{ir},'lineWidth',1);
            plot(t,pMedianLog-pRangeLog(2,:),'lineStyle',':','color',colours{ir},'lineWidth',1);
        end
        
    % 2. Lin-log scale  . . . . . . . . . . . . . . . . . . . . . . . . . .
    elseif ( strcmp(fig.yscale,'log') & strcmp(fig.xscale,'lin') )
        %ha{ir} = shadedErrorBar(t,pMedianLog, pRangeLog ,{['-',colours{ir}],'lineWidth',2},1);
        ha{ir} = plot(t,pMedianLog          ,'lineStyle','-','color',colours{ir},'lineWidth',2);
    
    % 3. Lin-lin scale  . . . . . . . . . . . . . . . . . . . . . . . . . .  
    elseif ( strcmp(fig.yscale,'lin') & strcmp(fig.xscale,'lin') )
        ha{ir} = plot(t,pMedian,'lineStyle','-','color',colours{ir},'lineWidth',2);
    else
        fprintf(1,'This yscale/xscale-combo is not an option...\n.')
    end
    
    % Annotate curves
    if fig.annotMedCurves
        xtxt     = 13;
        ytxt     = pMedianLog(end)+.15;
        %annotTxt = sprintf('n=%i',size(sAmps,1));
        %annotTxt = sprintf('%3.1f<=M<%3.1f, n=%i',mRanges(ir,:),size(sAmps,1));
        annotTxt = sprintf('%3.1f<=M<%3.1f',mRanges(ir,:));
        text(xtxt,ytxt,annotTxt,'fontSize',ftSize);
        %text(xtxt,ytxt,annotTxt,'fontSize',ftSize,'BackgroundColor','w');
    end
    
    % Save all values for easy reproduction
    outMat{ir}.t          = t;
    outMat{ir}.pMedianLog = pMedianLog;
    outMat{ir}.colours    = colours;
end 


%_ ADD NOISE PGD(t) CURVES  . . . . . . . . . . . . . . . . . . . . . . . .
if fig.plotNoiseAmps
    if     strcmp(band,'Pd'); noiseAmps=cellfun(@(x) x.pdN, trList.var5,'uniformOutput',0);
    elseif strcmp(band,'Pa'); noiseAmps=cellfun(@(x) x.paN, trList.var5,'uniformOutput',0);
    else   fprintf(1,'No noise amps available for chosen band.\n')
    end
    %     if     strcmp(band,'Pd'); noiseAmps = trList.var4;
    %     elseif strcmp(band,'Pa'); noiseAmps = trList.var5;
    %     else   fprintf(1,'No noise amps available for chosen band.\n')
    %     end
    
    % Make all signal traces have same length:
    % (1) if they are longer  than tRange(2), truncate them
    % (2) if they are shorter than tRange(2), replicate the last values
    nmax = floor(tRange(2)/snpLength);
    ntrn  = numel(noiseAmps);
    N    = zeros(ntrn,nmax);
    
    for itr = 1:ntrn
        ns = length(noiseAmps{itr});
        if ns>=nmax; N(itr,:) = noiseAmps{itr}(1:nmax);                                  % (1)
        else         N(itr,:) = [noiseAmps{itr}, repmat(noiseAmps{itr}(end),1,nmax-ns)];      % (2)
        end
    end
    tn = (1:length(N(1,:)))'*snpLength;
    
    % Compute percentiles
    pnl = prctile(N,pctl(1));
    pnm = prctile(N,pctl(2));
    pnu = prctile(N,pctl(3));
    
    hn = plot(tn,log10(pnm),'lineStyle','-','color',[.8 .8 .8],'lineWidth',2); 
end


%_ FORMAT FIGURE  . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
xlabel('Time since P-onset [sec]','fontSize',ftSize)
set(gca,'fontSize',ftSize,'xlim',[tmin tmax])
%if strcmp(fig.xscale,'log'); set(gca,'xscale','log','ylim',[-8 -1]); end
if strcmp(fig.xscale,'log') &&strcmp(band,'Pd'); set(gca,'xscale','log','ylim',[-8 -.5]); end
if strcmp(fig.xscale,'log') &&strcmp(band,'Pa'); set(gca,'xscale','log','ylim',[-4  1 ]); end
%if strcmp(fig.xscale,'log'); set(gca,'xscale','log','ylim',[-9 0]); end
if hasP; set(gca,'xscale','log','ylim',[-11 -5]); end

if strcmp(fig.yscale,'log')
    if strcmp(band,'Pd'); ylabel(['log_{10} (pgd) [m]']    ,'fontSize',ftSize); end
    if strcmp(band,'Pv'); ylabel(['log_{10} (pgv) [m/s]']  ,'fontSize',ftSize); end
    if strcmp(band,'Pa'); ylabel(['log_{10} (pga) [m/s^2]'],'fontSize',ftSize); end
else
    if strcmp(band,'Pd'); ylabel([band,' [m]'],'fontSize',ftSize); end
    if strcmp(band,'Pv'); ylabel([band,' [m]'],'fontSize',ftSize); end
    if strcmp(band,'Pa'); ylabel([band,' [m]'],'fontSize',ftSize); end
end
if fig.title
    title(sprintf('%s %s on records with %i<=R<%ikm, SNR>=%i, %s records, n_{pdc}=%3.1f',ornt,band,rRange,fig.snrMin,titleAppendix,fig.pdcExponent))
end

%fprintf(1,'\nTemporarily fixed ylim to [-7 -4]. Remove.\n\n')
%set(gca,'ylim',[-6 -2],'fontSize',ftSize)

if isnumeric(band);
    ylabel('log_{10} (PGV^{nb}) [m/s]','fontSize',ftSize)
    title(sprintf('PGV^{nb} at %3.1fHz - %3.1fHz on records with %i<=r<%ikm',fc(band,:),rRange))
end
posMain = get(gca,'position');
ylm     = get(gca,'ylim'); 
line([0 0],ylm,'color','k','lineWidth',1)



%_ PLOT THEORETICAL S-ARRIVALS    . . . . . . . . . . . . . . . . . . . . .
if ~strcmp(fig.xscale,'lin')
    vs  = 3.4;
    fprintf(1,['Plotting theoretical s-wave arrival assuming (hard-coded) vs = ',num2str(vs),'km/s...'])
    ylm = get(gca,'ylim');
    dy  = 0.1;
    for ir = 1:nr
        ts_arrival = trList.hypDist(idx{ir})/vs;
        plot(ts_arrival,ylm(1)+(nr-ir+1)*dy,'.','color',colours{ir},'markerSize',10)
        outMat{ir}.ts_arrival = ts_arrival;
    end
    fprintf(1,'done.\n')
end


%_ LEGEND    . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
if fig.legend
    posbak = get(gca,'position');
    
    lgdHandles = zeros(nr+1,1);
    lgdStrings = cell(nr+1,1);
    for ir = 1:nr
        if strcmp(fig.xscale,'log') |( strcmp(fig.xscale,'lin') &strcmp(fig.yscale,'lin') );
              lgdHandles(ir) = ha(ir);
        else  lgdHandles{ir} = ha{ir}.mainLine;
        end
        txt            = sprintf('%s: %3.1f <= M <%3.1f  (n=%i)',bt{ir},mRanges(ir,:),ntr(ir));
        lgdStrings{ir} = txt;
    end
    
    if fig.plotNoiseAmps; lgdHandles(end) = hn;
                          lgdStrings{end} = 'unscaled ambient noise';
    else                  lgdHandles      = lgdHandles(1:end-1);
                          lgdStrings      = lgdStrings(1:end-1);
    end
    
    % Legend (poition will be changed after hypocentral distances)
    hl   = legend(lgdHandles,lgdStrings);
    pos1 = get(gca,'position');
    
    %if strcmp(fig.xscale,'log') set(hl,'location','SouthEast','fontSize',ftSize)    
    if strcmp(fig.xscale,'log') set(hl,'location','SouthEastOutside','fontSize',ftSize)    
    else                        set(hl,'location','NorthWest'       ,'fontSize',ftSize)    
    end
    set(gca,'position',posbak)
end


%_ PLOT SLOPE SCALE
if fig.plotSlopeScale

    % Slope at bottom center of plot
    %     x0    = .07;
    %     xE    = .3;
    %     y0    = -8.7;
    %     xx    = x0:snpLength:xE;
    % 
    %     slopes = [0, 1.5, 3, 4, 5, 6];
    %     for is = 1:numel(slopes)
    %         c    = y0 - slopes(is)*log10(x0);
    %         logy = c  + slopes(is)*log10(xx);
    %         plot(xx,logy,'-k')
    %         text(xx(end)+.01,logy(end),num2str(slopes(is))  ,'fontSize',ftSize)
    %     end
    
    % Slope at top left of plot
    x0     = .015;
    xE     = .1;
    y0     = -1;
    xx     = x0:snpLength:xE;    
    slopes = [0, 1, 2, 3, 4, 5];
    
    for is = 1:numel(slopes)
        c    = y0 - slopes(is)*log10(xE);
        logy = c  + slopes(is)*log10(xx);
        plot(xx,logy,'-k')
        text(xx(1)-.004,logy(1),num2str(slopes(is))  ,'fontSize',ftSize)
    end
end






%% %%%%%%%%%%%%
%  4. KS TEST %
%%%%%%%%%%%%%%%
fprintf(1,'Computing and plotting KS tests ...')

[~,idxMin] = min(cellfun(@(x) size(x,2), amps));   % Find and choose shortest time series
t          = tvect{idxMin};                        % time vector of all bins
t          = t(t<=fig.tmaxKS);
nt         = numel(t);

%fprintf(1,'\nNOTE: KS-test is performed on log-amplitudes,\n      t-test  is performed on linear amplitudes\n')
%fprintf(1,sprintf('Using %4.2f as significance level\n',plevel))

% Compare distributions of consecutive pairs
ks   = cell (nr-1,1);
%st   = cell (nr-1,1);
tsep = zeros(nr-1,1);
isep = zeros(nr-1,1);
for ir = 1:nr-1
    
    fprintf(1,sprintf(' %i/%i ..',ir,nr-1))
    
    % Linear and logarithmic amlitudes of two consecutive bins
    if fig.considerPxDist
        A1_log = log10(pdKS{ir}(:,1:1000));
        A2_log = log10(pdKS{ir+1}(:,1:1000));
    else
        A1     = amps{ir};
        A2     = amps{ir+1};
        A1_log = log10(A1);
        A2_log = log10(A2); 
    end 
    
    
    % (1) KS test can be done on linear amplitudes
    % (2) T-test assumes normal distribution, i.e. must be performed on
    %     log-amplitudes
    pks = zeros(nt,1);
    prs = zeros(nt,1);
    %pst = zeros(nt,1);
    for it = 1:nt
    
        %         fprintf(1,'excluding extreme values from test\n')
        %         % Exclude outliers
        %         aa1  = A1(:,it);
        %         pct1 = prctile(aa1,[5 95]);
        %         aa11 = aa1(aa1>=pct1(1) & aa1<pct1(2));
        %         aa2  = A2(:,it);
        %         pct2 = prctile(aa2,[5 95]);
        %         aa22 = aa2(aa2>=pct2(1) & aa2<pct2(2));
        %         [~,pks(it)] = kstest2(aa11, aa22);
        
        [~,pks(it)] = kstest2(A1(:,it), A2(:,it));
        prs(it)     = ranksum(A1(:,it), A2(:,it));
        %[~,pks(it)] = kstest2(A1_log(:,it), A2_log(:,it));
        %[~,pst(it)] = ttest2 (A1_log(:,it), A2_log(:,it),'Vartype','unequal');
    end
    
    ks{ir} = pks; 
    rs{ir} = prs;
    %ks{ir} = pst;
    %st{ir} = pst;

    % Show distributions and percentiles for each time step
    %     for it = 1:nt
    %         
    %         figure(223); clf; 
    %         subplot(2,1,1); hold on; grid on; box on;
    %         plot(A1(:,it),1  ,'.k')
    %         plot(A2(:,it),1.1,'.r')
    %         prct1 = prctile(A1(:,it),[0,1,5,16,50,84,95,99,100]);
    %         prct2 = prctile(A2(:,it),[0,1,5,16,50,84,95,99,100]);
    %         plot(prct1,1  ,'vk','markerSize',12)
    %         plot(prct2,1.1,'vr','markerSize',12)
    %         %set(gca,'xlim',[1e-8 1e-1],'xscale','log','ylim',[0.5 1.5])
    %         set(gca,'xscale','lin','ylim',[0.5 1.5])
    %     
    %         subplot(2,1,2); hold on; grid on; box on;
    %         p = pks(it);
    %         if p<1e-3; p=1e-3; end
    %         plot(1,p,'ok','markerFaceColor','r','markerSize',10)
    %         set(gca,'xlim',[0.9 1.1],'ylim',[1e-3 1],'yscale','log')
    %         1+1;
    %     end

    %     % Subsample larger array so that they have same size
    %     nit = 100;
    %     nA1 = size(A1,1);
    %     nA2 = size(A2,1);
    %     pksMat = zeros(nt,nit);
    %     A1master = A1;
    %     A2master = A2;
    %     for iit = 1:nit
    %         if nA1>nA2; irand = randi(nA1,nA2,1); A1 = A1master(irand,:); end
    %         if nA2>nA1; irand = randi(nA2,nA1,1); A2 = A2master(irand,:); end
    %         pks = zeros(nt,1);
    %         for it = 1:nt
    %             [~,pks(it)] = kstest2(A1(:,it), A2(:,it));
    %         end
    %         pksMat(:,iit)=pks;
    %     end
    %     ks{ir} = median(pksMat,2);

    % tsep: compute time at which <pks> goes (and stays) below <plevel> 
    %idx_sep = find(pst>=plevel,1,'last')+1;
    %idx_sep = find(pks>=plevel,1,'last')+1;
    idx_sep = find(prs>=plevel,1,'last')+1;

    if isempty(idx_sep); idx_sep=1; end
    if idx_sep>nt;       idx_sep=nt; end
    
    tsep(ir) = t(idx_sep);
    isep(ir) = idx_sep;
end

% Are there enough samples for KS test?
%n3               = numel(a3);
%n4               = numel(a4);
%enough_samples_q = (n3*n4)/(n3 + n4);
%txt = sprintf('%4.1f > 4?',enough_samples_q);
%text(.8*tmax,.8,txt,'fontSize',ftSize,'BackgroundColor','w');


% Plot tests
if ~fig.multipleKSplots; tAx = subplot(nlines,1,5); hold on; box on; grid on; 
else                     tAx = zeros(nr-1,1);
end

for ir = nr-1:-1:1

    if fig.multipleKSplots; subplot(nlines,1,4+ir); hold on; box on; grid on; end
    
    c1 = colours{ir};
    c2 = colours{ir+1};
    
    hks   = plot(t',rs{ir},'-','color',c1,'lineWidth',fig.ksLineWidth);
            plot(t',rs{ir},'-','color',c2,'lineWidth',fig.ksLineWidth,'lineStyle','-.');
    htsep = plot(tsep(ir),rs{ir}(isep(ir)),'ow','lineWidth',2);
            plot(tsep(ir),rs{ir}(isep(ir)),'ok','lineWidth',1);
    %hst = plot(t',st{ir},'-b'                 ,'lineWidth',fig.ksLineWidth);
    set(gca,'xlim',[tmin tmax],'ylim',[0,1],'fontSize',ftSize);
    ylabel('KS p-value','fontSize',ftSize)

    % Plot duration CDFS of median M in lower M-bin
    %pd = plot(tcdfDuration{ir+1},ecdfDuration{ir+1},'--r');
    %durationText = sprintf('M_{med}%3.1f',M_cdf(ir+1));
    %xtxt = tcdfDuration{ir+1}(round(nsample/2))+.8;
    %ytxt = 0.4;
    %text(xtxt,ytxt,.7,durationText,'fontSize',ftSize,'color','r','BackgroundColor','w');

    % Plot tsep
    %plot(tsep(ir), ks{ir}(isep(ir)),'or')

    % Annotate figures
    annotText = sprintf('%s vs   %s',bt{ir+1}, bt{ir});
    %     if ( rem(mRanges(ir,1),1)==0 && rem(mRanges(ir+1,1),1)==0 ) ; % then its an integer
    %     else      annotText = sprintf('M%3.1f+ vs M%3.1f+',mRanges(ir+1,1),mRanges(ir,1));
    %     end
    
    if strcmp(fig.xscale,'lin'); xtxt = .6*tmax; end
    if strcmp(fig.xscale,'log'); xtxt = .1*tmax; end
    ytxt = .8;
    %text(xtxt,ytxt,annotText,'fontSize',ftSize,'BackgroundColor','w');
    ylm = get(gca,'ylim');
    line([0 0],ylm,'color','k','lineWidth',1)

    %if ir~=nr-1; set(gca,'xTickLabel',''); end
    if ir==1 &&fig.multipleKSplots; pos1 = get(gca,'position'); ht = pos1(2)+pos1(4); end
    if ir==nr-1; pos2 = get(gca,'position'); hb = pos2(2);         end
    if strcmp(fig.xscale,'log'); set(gca,'xscale','log'); end
    pos = get(gca,'position'); pos(3)=posbak(3); set(gca,'position',pos)
end
xlabel('Time since P-onset [sec]','fontSize',ftSize)
fprintf(1,'done.\n')
xlm = get(gca,'xlim');

if isfield(fig,'yscale_KS')
    if strcmp(fig.yscale_KS,'log')
    ytx = logspace(-3,0,4);
    set(gca,'yscale',fig.yscale_KS,'ylim',[1e-3, 1],'ytick',ytx)
    end
    pos    = get(gca,'position');
    pos(2) = .85*pos(2);
    set(gca,'position',pos)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%  5. Plot Duration Scalings %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fig.plotDurationScalings
    mm = [3, 4, 5, 6, 7, 8];
    nm = numel(mm);
    n  = 1e4;    % No. of random samples for computing durations
    
    % Rupture duration [sec] & length [km]
    [D,L,c] = get_eqSourceDuration_WellsCoppersmith94(mm,n);
    
    % Subplot
    if fig.multipleKSplots; lineNo      = 4+nr;
    else                    lineNo      = 6;
    end
    subplot(nlines,1,lineNo); hold on; box on; grid on;
    
    for im = 1:nm
        
        % Plot DURATION cdf
        [fd,xd] = ecdf(D(:,im));
        plot(xd,fd, 'lineWidth',1,'color','k')
        
        % Annotate
        txt = sprintf('M%3.1f',mm(im));
        xt  = median(xd)-.4*median(xd);
        yt  = .82;
        text(xt,yt,txt,'fontSize',ftSize,'BackgroundColor','w')
    end
    set(gca,'xscale','log', 'fontSize',ftSize)
    
    % title(sprintf('Source scalings of Wells and Coppersmith, 1994, BSSA, using a=%4.2f & b=%4.2f, vr=unif[2.4 3.0]km/s ', ...
    %     c.a, c.b), 'fontSize',ftSize)
    xlabel(gca,'Rupture Duration Scaling [sec]', 'fontSize',ftSize)
    ylabel(gca,'p(D''<=D)', 'fontSize',ftSize)
    pos = get(gca,'position'); pos(3)=posbak(3); set(gca,'position',pos)
    set(gca,'xlim',xlm)
end

% Modify positions
% h   = ht-hb;                        % Combined height of all test-subplots
% nh  = 4;                            % Subplot has nh times the height of space between subplots
% dys = h / ((nr-1)*nh + (nr-2));     % Height of space between subplots
% dyl = nh*dys;                       % Heigth of subplots
% for ir = 1:nr-1
%     subplot(tAx(ir));
%     pos    = get(gca,'position');
%     ypos   = ht-ir*dyl-(ir-1)*dys;
%     pos(2) = ypos;
%     pos(4) = dyl;
%     set(gca,'position',pos)
%     if ir==1;    yposTopFirst = ypos+dyl; end   % Save it for shifting legend
%     if ir==nr-1; yposLast     = ypos;     end   % Save it for shifting distance-cdf
% end





%% %%%%%%%%%%%%%%%%%%%%%%%%%
%  5. Plot observed slopes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fig.plotObsSlopes
   
    hsp(2)      = subplot(2,1,2); hold on; grid on; box on;
    pos1        = get(hsp(1),'position');
    pos2        = get(hsp(2),'position');
    pos2([1,3]) = pos1([1,3]);
    set(gca,'position',pos2)
    
    
    % Find minimum common length of tdpa-curves in all mBins
    nsmin = zeros(nr,1);
    for ir = 1:nr
        nsmin(ir) = size(signal{ir}.amps,2);
    end
    
    nsmin = min(nsmin);
    slopeMat = zeros(nr,nsmin-1);
    
    for ir = 1:nr
        
        % Unpack signal data
        ts  = signal{ir}.t';
        psm = prctile(signal{ir}.amps,pctl(2));
        
        
        % Compute slope on log-log plot
        tlog  = log10(ts);
        alog  = log10(psm);
        da    = diff(alog);
        dt    = diff(tlog);
        slope = da./dt;
        
        slopeMat(ir,:) = slope;
        slopeSmooth = movingAverage(slope,0);
        
        subplot(hsp(2)); hold on;
        plot(tlog(1:end-1),slopeSmooth,'r','color',trList{ir},'lineWidth',1)
        set(gca,'ylim',[0 9])
    end
    set(gca,'fontSize',ftSize,'xlim',[log10(tmin) log10(tmax)])
end



%% %%%%%%%%%%%%%%%%%
%  6. Print Figure %
%%%%%%%%%%%%%%%%%%%%

if fig.printPng |fig.printEps
    
    fprintf(1,'Printing figures ..')
    set(gcf,'PaperPositionMode','auto')
    
    % Automatic naming commands moved to APPENDIX
    
    
    if fig.printPng; print('-dpng'  ,[fig.printFullName,'.png']); end
    if fig.printEps; print('-depsc2',[fig.printFullName,'.eps']); end
    %print('-dpng','-r600',[outFullName,'.png'])
    %print('-dpdf',[outFullName,'.pdf'])
    fprintf(1,'done.\n')
end




%% %%%%%%%%%%%%%%%%%%%%%%%
%  7. Save Output Matrix %
%%%%%%%%%%%%%%%%%%%%%%%%%%
if fig.saveTdpaMat
    fprintf(1,'Saving tdpa-matrix ..')
    save(fig.tdpaMatFullName,'outMat')
    fprintf(1,'done.\n')
end



%% %%%%%%%%%%%%%
%  A. Appendix %
%%%%%%%%%%%%%%%%

%     nt  = size(sAmps,2);
%     psl = zeros(1,nt);
%     psm = zeros(1,nt);
%     psu = zeros(1,nt);
%     for it = 1:nt
%         
%         tmp = sAmps(:,it);
%         lgc = tmp~=0;
%         
%         % Only plot function for times when at least 10 non-zero
%         % measurements are available
%         if sum(lgc)>=10
%             psl(it) = prctile(tmp(lgc),pctl(1));
%             psm(it) = prctile(tmp(lgc),pctl(2));
%             psu(it) = prctile(tmp(lgc),pctl(3));
%         else
%             psl(it) = 0;
%             psm(it) = 0;
%             psu(it) = 0;
%         end
%     end

%outPath = sprintf('~/programs/filterBank/fig/i%i/cf_tdpa/new/',iN);
%mString = strrep(sprintf('%3.1f_',fliplr(mRanges(:,1)')),'.','p');
%mString = sprintf('%i',fliplr(mRanges(:,1)'));
%tString = sprintf('%i',tRange(2));
%     if ~isnumeric(band)
%         %outFullName = sprintf('~/programs/filterBank/fig/i%i/cf_tdpa/new/%s_%s_%i_to_%ikm_%srec',iN,band,ornt,rRange,titleAppendix);
%         outFullName = sprintf('~/programs/filterBank/fig/i%i/cf_tdpa/new/%s_%s_%i_to_%ikm_m%s_t%ssec_%srec_noTests', ...
%             iN,band,ornt,rRange,mString,tString,titleAppendix);
%     else
%         outFullName = sprintf('~/programs/filterBank/fig/i%i/cf_tdpa/new/band_%i_%s_%i_to_%ikm_%srec_noTests',iN,band,ornt,rRange,titleAppendix);
%     end
%
%     if strcmp(fig.xscale,'log'); outFullName = sprintf('%s_logt',outFullName); end
%     if strcmp(fig.xscale,'lin'); outFullName = sprintf('%s_lint',outFullName); end
