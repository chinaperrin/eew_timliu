o.useRconstraint = 0;
o.plotSSR=1;

testOutDirFullName = '~/programs/gba/tests/test3/out/';

for it = 1:ntss
    
    snippet = snpList(it);
    if ~o.quiet
        if (it==1); fprintf(1,['    ',num2str(ntss),' snippets: ',num2str(snippet),' .. '])
        else          fprintf(1,[num2str(snippet),' .. '])
        end
    end
    
    train.az = zamps{snippet};
    train.ah = hamps{snippet};
    [mle]    = estimate_params_SSR_i36(zTarget,hTarget,train,wti0,snippet,1:9,{'Z'},nsim,[],[]);
    
    % For each station
    outName     = sprintf('%s_snp%i',zTarget.stationName{1},snippet);
    outFullName = sprintf('%s%s',testOutDirFullName,outName);
    eval([outName,'.mle     = mle;']);
    eval([outName,'.zTarget = zTarget;']);
    %eval([outName,'.train   = train;']);
    
    save(outFullName,outName)
    
    % Save figure 
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2',[outFullName,'.eps'])
end


% Meta info
outDirFullName       = '/Users/mameier/programs/gba/tests/test3/out/';
test3OutFullName     = sprintf('%sgba_test3_metaInfo',outDirFullName);
gba_test3.zamps      = zamps;
gba_test3.hamps      = hamps;
gba_test3.mVect      = train.m;
gba_test3.rVect      = train.r;
gba_test3.mm         = mm;
gba_test3.rr         = rr;
gba_test3.snpLength  = snpLength;
gba_test3.fc         = fc;
gba_test3.fMode      = fMode;
gba_test3.fOrder     = fOrder;
save(test3OutFullName,'gba_test3')



%% APPENDIX
%dlmwrite('lmmn.txt',mle.lmmn)
%dlmwrite('lmrn.txt',mle.lmrn)    