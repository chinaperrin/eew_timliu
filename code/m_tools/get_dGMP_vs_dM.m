function [hf,pct] = get_DGMP_vs_dM(trList,methPrime,gmMeasure,twPrime,opts)

global ftSize
ntr = numel(trList.m);
np  = numel(opts.compPctl);

% Unpack opts
tms       = opts.tms;
twMaster  = opts.twMaster;
crange    = opts.crange;
clrNUM    = opts.clrNUM;
clrMAP    = opts.clrMAP;
twPlotMin = opts.tRange(1);
twPlotMax = opts.tRange(2);
gmMin     = opts.yRange(1);
gmMax     = opts.yRange(2);

estimatorName = 'hat';
methName      = strrep(methPrime,'_ss','_{ss}');
lnWidth_ss    = 1;
dc            = 1;  % Default color increment

% Find time-index of entry that corresponds to warning time of <twPrime> 
[val,itw] = min(abs(twMaster-twPrime));

% Observed and predicted ground motion
if     strcmp(gmMeasure,'immi'); gmObs = cellfun(@(x) x.zen, trList.intensity);
elseif strcmp(gmMeasure,'pga') ; gmObs = trList.pga;
elseif strcmp(gmMeasure,'pgv') ; gmObs = trList.pgv;
elseif strcmp(gmMeasure,'m')   ; gmObs = trList.m;
elseif strcmp(gmMeasure,'rh')  ; gmObs = trList.hypDist;
                                 dc    = 10;
else   error('This GM measure has not yet been coded up.')
end

% Sort list: Make traces with largest observed GM come last so that they 
% will be plotted on top
if ~strcmp(gmMeasure,'rh'); [gmObs,idx] = sort(gmObs);
else                        [gmObs,idx] = sort(gmObs,'descend');
end
trList.sortList(idx)

% Interpolate all residual curves onto same time vector
ntm  = numel(twMaster);
GMP  = zeros(ntr,ntm);
dGMP = zeros(ntr,ntm);
MHAT = zeros(ntr,ntm);
dM   = zeros(ntr,ntm);

hf=figure(767); clf; whitebg('w'); hold on; box on; grid on;
hasEstimate=true(ntr,1);
for itr = 1:ntr
    
    ssrte = trList.var1{itr};
    mObs  = trList.m(itr);
%     if isempty(regexp(methPrime,'pythia','once')); mObs  = trList.m(itr);
%     else                                           mObs  = 0;
%     end
    
    % Warning time
    te   = ssrte.gmpp.te; 
    twms = -1*(tms-te);
    
    % Extract intensity prediction  . . . . . . . . . . . . . . . . . . . .
    eval(sprintf('gmp  = ssrte.%s.%s.%s;',methPrime,gmMeasure,estimatorName))

    if isempty(regexp(methPrime,'pythia','once')); 
        eval(sprintf('mHat = ssrte.%s.%s.%s;',methPrime,'m' ,estimatorName))
    else    
        mHat = mObs*ones(size(tms));
    end

    if ~isempty(gmp)
        %if ~isempty(regexp('PGA','gmMeasure','once')); gmp(gmp~=0)=10.^gmp(gmp~=0); end
        
        gmp_intpol                = interp1(twms,gmp,twMaster);        % Interpolated to common time vector
        iLastVal                  = find(~isnan(gmp_intpol),1,'last'); % Replace values after last update with last predictions
        gmp_intpol(iLastVal:end)  = gmp_intpol(iLastVal);              % Replace zeros from before first predictions with nans
        gmp_intpol(gmp_intpol==0) = nan;                               % so that they will not be plotted
        
        mhat_intpol                 = interp1(twms,mHat,twMaster);        % Interpolated to common time vector
        iLastVal                    = find(~isnan(mhat_intpol),1,'last'); % Replace values after last update with last predictions
        mhat_intpol(iLastVal:end)   = mhat_intpol(iLastVal);              % Replace zeros from before first predictions with nans
        mhat_intpol(mhat_intpol==0) = nan;                               % so that they will not be plotted
        
        % GMP & GMP error
        GMP (itr,:) = gmp_intpol;
        MHAT(itr,:) = mhat_intpol;
        dM  (itr,:) = mObs-mhat_intpol;
        if isempty(regexp('PGA','gmMeasure','once')); dGMP(itr,:) = gmObs(itr)-gmp_intpol;
        else                                          dGMP(itr,:) = log10(gmObs(itr)./gmp_intpol);
                                                      error('not yet adapted to use any other gmMeasure than IMMI\n')
        end
        
        if isempty(regexp('PGA',gmMeasure,'once')); clrID = round(clrNUM*(      gmObs(itr) -crange(1))/range(crange));
        else                                        clrID = round(clrNUM*(log10(gmObs(itr))-crange(1))/range(crange));
        end
        if clrID<1;      clrID=1;      end
        if clrID>clrNUM; clrID=clrNUM; end
        lnColor = clrMAP(clrID,:);
        
        plot(dM(itr,itw),dGMP(itr,itw),'ok','markerFaceColor',lnColor,'lineWidth',lnWidth_ss)
    else
        hasEstimate(itr)=false;
    end
end
fprintf(1,'%i/%i have no estimate.\n',ntr-sum(hasEstimate),ntr)
GMP  = GMP (hasEstimate,:);
dGMP = dGMP(hasEstimate,:);
ntr  = size(GMP,1);


% % Compute percentiles
% pct = zeros(np,ntm);
% for ipt = 1:np
%     pct(ipt,:) = prctile(dGMP,opts.compPctl(ipt));
%     if opts.compPctl(ipt)==50; lnWidth=2; else lnWidth=1; end
%     if opts.plotPctl(ipt); plot(twMaster,pct(ipt,:),'-k','lineWidth',lnWidth); end
% end
fprintf(1,'function not yet adapted to compute percentiles\n')
pct=[];

%set(gca,'ylim',[gmMin gmMax],'xlim',[twPlotMin,twPlotMax],'fontSize',ftSize,'FontName','Avenir','xdir','reverse');
set(gca,'fontSize',ftSize,'FontName','Avenir');
line(get(gca,'xlim'), [0, 0],'lineWidth',1,'color','k','lineStyle','-.')
%line([0 0]          , get(gca,'ylim')       ,'lineWidth',1,'color','k','lineStyle','-' )
%title(sprintf('Sites exceeding MMI%i',immiPrime),'fontSize',ftSize,'FontName','Avenir')
xlabel('Magnitude Residual','fontSize',ftSize,'FontName','Avenir')
ylabel('GMP Residual','fontSize',ftSize,'FontName','Avenir')
title(sprintf('dGMP vs dM at %is before strong GM onset; from method %s',twPrime,methName),'fontSize',ftSize,'FontName','Avenir');

colormap(clrMAP)
hc1 = colorbar;
ctx = crange(1):dc:crange(end);
xtx = (ctx-crange(1))./range(crange);
nt  = numel(xtx);
set(hc1,'xtick',xtx);

if strcmp(gmMeasure,'immi')
    hc1.Label.String='MMI_{obs}';
    ylabel('MMI_{obs}-MMI_{pred}','fontSize',ftSize,'FontName','Avenir')
    xtL = {linspace(crange(1),crange(end),nt)};
elseif strcmp(gmMeasure,'pga')
    %set(gca,'ylim',[gmMin gmMax])
    hc1.Label.String='PGA_{obs} [m/s^2]';
    ylabel('log_{10}[PGA_{obs}/PGA_{pred}]','fontSize',ftSize,'FontName','Avenir')
    xtL = {10.^ctx};
elseif strcmp(gmMeasure,'pgv')
    hc1.Label.String='PGV_{obs} [m/s]';
    ylabel('log_{10}[PGV_{obs}/PGV_{pred}]','fontSize',ftSize,'FontName','Avenir')
    xtL = {10.^ctx};
elseif strcmp(gmMeasure,'m')
    hc1.Label.String='M_{catalog}';
    ylabel('M_{catalog}-M_{est}','fontSize',ftSize,'FontName','Avenir')
    xtL = {linspace(crange(1),crange(end),nt)};
elseif strcmp(gmMeasure,'rh')
    hc1.Label.String='Rh_{catalog}';
    ylabel('Rh_{catalog}-Rh_{est}','fontSize',ftSize,'FontName','Avenir')
    xtL = {linspace(crange(1),crange(end),nt)};
end
set(hc1,'xTickLabel',xtL,'fontSize',ftSize,'FontName','Avenir');

if opts.printPng
    figFullName = sprintf('%s_%s_%s',opts.figFullName,methPrime,gmMeasure);
    set(gcf,'PaperPositionMode','auto');
    print('-dpng',figFullName)
end
