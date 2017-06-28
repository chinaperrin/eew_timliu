% Save random number generator settings
    %seed_rng_stacks = rng;
    %seedFullName    = sprintf('~/programs/filterBank/var/wform_stacks/rng_seed_i%i_%intr.mat',iN,numel(trList.m));
    %save(seedFullName,'seed_rng_stacks')
    rng(seed_rng_stacks) % --> reset rng settings
    %rand(1,10,1)'
    
    %mRanges = [6.5 7; 6 6.5; 5.5 6; 5 5.5; 4.5 5];
    mRanges = [6.5 7.7; 6 6.5; 5.5 6; 5 5.5; 4.5 5; 4 4.5; 3.5 4; 3 3.5; 2.5 3];
    rRange  = [0 25];
    
    % Exclude SM records for M<5
    idx_M5sm      = find(zList.m<5 & ~strcmp(zList.instrCode,'H')); 
    nz            = numel(zList.m);
    idx_noSmallSM = setdiff(1:nz,idx_M5sm);
    allList       = zList.selectSubList(idx_noSmallSM);
    %allList       = hList.selectSubList(idx_noSmallSM);

    % Exclude NGA records
    %nga_idx    = find(strcmp(allList.dataSetName,'ngawest1'));
    %nonga_idx  = setdiff(1:numel(allList.m),nga_idx);
    %allList    = allList.selectSubList(nonga_idx);

    
    trList  = allList.selectSubList( find(allList.hypDist>=rRange(1) & allList.hypDist<rRange(2) & allList.snr>50) );
    
    opts.process        = 1;
    opts.intMode        = 'afterPx';
    opts.verbose        = false;
    opts.scp_wforms     = true;
    opts.saveStacks     = true;
    opts.plotNormStacks = true;
    opts.parallel       = true;
    opts.filterStacks   = false;
    opts.fUpStack       = 10;
    
    % Restore rand-no.-generator settings
    rng(s);

    ntrmax    = 1000;
    ns_prePx  = 100;
    ns_postPx = 200;
    ns_stack  = ns_prePx + 1 + ns_postPx;
    idx_stack = -ns_prePx:1:ns_postPx;
    
    nm        = size(mRanges,1);
    allStacks = cell(nm,1);
    allCounts = cell(nm,1);
    allNtr    = zeros(nm,1);
    
    % Parallel
    if opts.parallel
        uname = getenv('USER');
        if strcmp(uname,'mameier'); nWorkers = 2; else nWorkers = 24; end
        if matlabpool('size')==0; matlabpool('open',nWorkers);          % checking to see if my pool is already open
        else                      fprintf(1,['Matlabpool already open (',num2str(matlabpool('size')),' workers)\n'])
        end
    end
    
    
    

    for im = 1:nm
        
        mLo     = mRanges(im,1);
        mUp     = mRanges(im,2);
        tmpList = trList.selectSubList( find(trList.m>=mLo &trList.m<mUp) );
        ntr     = numel(tmpList.m);
        fprintf(1,sprintf('\n\nM: %3.1f - %3.1f\n',mLo,mUp))
        
        % Shorten list if it has more than ntrmax entries
        if ntr>ntrmax
            idxRand = randsample(ntr,ntrmax);
            tmpList = tmpList.selectSubList(idxRand);
            ntr     = numel(tmpList.m);
        end
        
        [allStacks{im},allCounts{im}] = stack_wforms_around_px(tmpList,ns_prePx,ns_postPx,opts);
        allNtr(im)                    = ntr;
    end
    
    if opts.saveStacks
        stackOut.allStacks = allStacks;
        stackOut.allCounts = allCounts;
        stackOut.allNtr    = allNtr;
        stackOutFullName   = sprintf('~/programs/filterBank/var/wform_stacks/i%i_%intr.mat',iN,numel(trList.m));
        save(stackOutFullName,'stackOut')
    end
    
    if opts.plotNormStacks
       
        % Unpack saved info if script was run on bigstar
        %load ~/programs/filterBank/var/wform_stacks/i35_6821ntr.mat
        %allStacks = stackOut.allStacks;
        %allCounts = stackOut.allCounts;
        %allNtr    = stackOut.allNtr;
        
        figure(1); clf; hold on; grid on; box on;
        %colours = {'r','k','m','c','g','y','b','r','k','m','c','g','y','b','r'};
        colours = {'r','k','m','c','g','y','b',[0 .4 0],[.5 0 0]};
        %lnStyle = {'-r','-k','-m','-c','-g','-y','-b',':r',':k',':m',':c',':g',':y',':b',':r'};
        %set(gca,'xlim',[-50 50],'ylim',[0 .1])
        
        normStacksAcc = zeros(nm,ns_stack);
        normStacksVel = zeros(nm,ns_stack);
        normStacksDsp = zeros(nm,ns_stack);
        for im = 1:nm
            
            mUp = mRanges(im,2);
            if mUp>=5; lnStyle = '-'; 
            else       lnStyle = ':';
            end
            
            % ACCELERATION
            % Normalise acc-stacks to [0 1]
            iStack              = allStacks{im}.acc;
            yr                  = max(iStack) - min(iStack);
            normStacksAcc(im,:) = (iStack-min(iStack))./yr;
            
            if opts.filterStacks;   fprintf(1,'Lowpass-filtering stack\n')
                                    normStacksAcc(im,:) = bworth(normStacksAcc(im,:),100,opts.fUpStack,'low',4,'causal');
            else                    fprintf(1,'Not applying any filter to stack\n')
            end
            
            %plot(idx_stack,normStacksAcc(im,:),'lineStyle',lnStyle,'color',colours{im},'lineWidth',2)
            plot(idx_stack,iStack,'lineStyle',lnStyle,'color',colours{im},'lineWidth',2)
            ylm = get(gca,'ylim');
            line([0 0],ylm,'color','r','lineWidth',2)
            
        end
    end
    
    if opts.parallel; matlabpool close; end
        