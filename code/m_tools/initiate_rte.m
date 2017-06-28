function rte = initiate_rte(methodList,ntargets,ntimesteps,nit,opts)

% Initiate real-time estimation structure for all methods in methodList
% Structure: rte.methodName.quantity.fields

% ntimesteps = number of time steps
% ntargets   = number of target stations

fprintf(1,'Note: single-site methods need to have ''_ss'' in method name\n')
nm            = numel(methodList);
flg_plotHat   = opts.plotHat;
flg_plotRange = opts.plotRange;
%if isempty(rtePlotFlags); rtePlotFlags=false(nm,1); end

% Some quantities have one estimation per earthquake, e.g. magnitude estimate
eqQuantity.hat   = zeros(1,ntimesteps,'single');
eqQuantity.hat2  = zeros(1,ntimesteps,'single');
eqQuantity.mean  = zeros(1,ntimesteps,'single');
eqQuantity.var   = zeros(1,ntimesteps,'single');
eqQuantity.vals  = cell (1,ntimesteps);
eqQuantity.prct  = cell (1,ntimesteps);
eqQuantity.llh   = cell (1,ntimesteps);

% Others are per target site, e.g. ground motion predictions
siteQuantity.hat   = zeros(ntargets,ntimesteps,'single');
siteQuantity.hat2  = zeros(ntargets,ntimesteps,'single');
siteQuantity.mean  = zeros(ntargets,ntimesteps,'single');
siteQuantity.var   = zeros(ntargets,ntimesteps,'single');
siteQuantity.vals  = cell (ntargets,ntimesteps);
siteQuantity.prct  = cell (ntargets,ntimesteps);
siteQuantity.llh   = cell (ntargets,ntimesteps);
%siteQuantity.pEx   = cell (ntargets,ntimesteps);
siteQuantity.pEx   = repmat({zeros(nit,ntimesteps,'single')},ntargets,1);

meta.longName = '';

opts.eval      = true;
opts.plotGMP   = true;
opts.plotSkill = false;
opts.plotRange = false;
opts.plotMean  = false;

for im = 1:nm
    mName = methodList{im};
    
    isSingleSiteMethod = logical(regexp(mName,'_ss')>0);

    % For single-station methods, all quantities are siteQuantitites, e.g.
    % magntudes estimates are computed for every
    if isSingleSiteMethod
        
        eval(sprintf('rte.%s.m    = siteQuantity;',mName));
        eval(sprintf('rte.%s.lat  = siteQuantity;',mName));
        eval(sprintf('rte.%s.lon  = siteQuantity;',mName));
        eval(sprintf('rte.%s.pga  = siteQuantity;',mName));
        eval(sprintf('rte.%s.pgv  = siteQuantity;',mName));
        eval(sprintf('rte.%s.immi = siteQuantity;',mName));
        eval(sprintf('rte.%s.rh   = siteQuantity;',mName));
        eval(sprintf('rte.%s.re   = siteQuantity;',mName));
    else
        % For multi-station methods, some quantities are only computed per
        % earthquake, not per site
        eval(sprintf('rte.%s.m    = eqQuantity;'  ,mName));
        eval(sprintf('rte.%s.lat  = eqQuantity;'  ,mName));
        eval(sprintf('rte.%s.lon  = eqQuantity;'  ,mName));
        eval(sprintf('rte.%s.pga  = siteQuantity;',mName));
        eval(sprintf('rte.%s.pgv  = siteQuantity;',mName));
        eval(sprintf('rte.%s.immi = siteQuantity;',mName));
        eval(sprintf('rte.%s.rh   = siteQuantity;',mName));
        eval(sprintf('rte.%s.re   = siteQuantity;',mName));
    end
    
    eval(sprintf('rte.%s.meta           = meta;'             ,mName));
    eval(sprintf('rte.%s.opts.eval      = true;'             ,mName));
    eval(sprintf('rte.%s.opts.plotHat   = flg_plotHat  (im);',mName));
    eval(sprintf('rte.%s.opts.plotRange = flg_plotRange(im);',mName));
end