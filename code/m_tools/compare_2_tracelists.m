dataSetNames = {'/scratch/memeier/data/ngawest/out/i34/causal_0p05/'; ...
                '/scratch/memeier/data/ngawest/out/i34b/causal_0p05/'};

% i34  = with    bworth pre-filtering       --> aList
% i34b = without bworth pre-filtering       --> bList
            
[TraceList,fc,snpLength] = import_trLists(dataSetNames,'traceList');     

fprintf(1,'\nSplitting TraceList into z-, e- & nList ... ')
idx_Z = find(strcmp(TraceList.orntCode,'Z'));
idx_E = find(strcmp(TraceList.orntCode,'E'));
idx_N = find(strcmp(TraceList.orntCode,'N'));

zList = TraceList.selectSubList(idx_Z);
eList = TraceList.selectSubList(idx_E);
nList = TraceList.selectSubList(idx_N);  

% Manually find index of when second list starts...
plot(datenum(zList.eqDate))

aList = zList.selectSubList(1:1409);
bList = zList.selectSubList(1410:2819);

% Use only traces with high SNR
aIdx  = find(aList.snr>100);
bIdx  = find(bList.snr>100);
aList = aList.selectSubList(aIdx);
bList = bList.selectSubList(bIdx);



%% Plot percentiles
for iband = 1:9
    aAmps = cell2mat(cellfun(@(x) x(iband,1:40)', aList.amax','UniformOutput',0));
    aPct = prctile(aAmps',[5, 50, 95]);
    
    bAmps = cell2mat(cellfun(@(x) x(iband,1:40)', bList.amax','UniformOutput',0));
    bPct = prctile(bAmps',[5, 50, 95]);
    
    figure(iband); clf; hold on; grid on
    p1=plot(aPct','r','lineWidth',2);
    p2=plot(bPct','k','lineWidth',2);
    set(gca,'yscale','log','fontSize',15);
    title(sprintf('%i^{th} band',iband),'fontSize',15)
    xlabel('Snippet','fontSize',15)
    ylabel('Peak amplitude','fontSize',15)
   
    if iband==1
        l1 = legend([p1(1);p2(1)],'with pre-filtering (5, 50, 95^{th} prctile)','wo/ pre-filtering (5, 50, 95^{th} prctile)');
        set(l1,'fontSize',12,'location','south')
    end    
    1+1;
end






%% Compare individual waveforms
na = numel(aList.m);

fLow_prefilt = 0.075;
fOrder       = 4;
fMode        = 'causal';
ntap         = 100;

for ia = 1:na
    
    % Find corresponding bList entry
    aFullName = aList.fullName{ia};
    ib        = find(cellfun(@(x) ~isempty(x), regexp(bList.fullName,aFullName)));
    
    ppxIdx    = aList.ppxIdx(ia);
    
    if ~isempty(ib)
        
        % Import aList-entry and prefilter
        [sraw,meta] = read_ascii_wform_nga(aFullName,1);
        ns          = numel(sraw);
        s           = sraw*9.80665; % From [g] to [ms^-2] (BIPM: 1g = 9.80665 ms^-2)    
        sr          = meta.sr;

        % a-type processing
        sm   = s - mean(s(1:ppxIdx));
        stap = taper(sm,sr,ntap);
        a    = bworth(stap,sr,fLow_prefilt,'high',fOrder,fMode);
        
        % b-type processing
        b  = taper(s,sr,ntap);
        
        % Plot both 
        figure(44); clf; hold on; grid on
        plot(a,'r','lineWidth',2)
        plot(b,'k')
        set(gca,'xlim',[1 500])
    
    else
        figure(44); clf; hold on; grid on
        fprintf(1,'no corresponding bList entry\n')
    end
    1+1;

end
