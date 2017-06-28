function [hf] = get_tsep_matrix(trList,band,mRanges,rRange,plevel,fig)

global fc ftSize iN snpLength

set(0,'DefaultFigureWindowStyle','docked')

nr      = size(mRanges,1);
hf      = figure;
nsample = 1e3;

if strcmp(trList.orntCode{1},'Z'); ornt = 'Vertical';
else                               ornt = 'Horizonal';
end 
 
testLineWidth = 2;
 
tRange = [0 100];
% tmin = tRange(1);
% tmax = tRange(2);
% if tmin==0 & strcmp(fig.xscale,'log'); tmin = 0.04; end

% Read peak amps from traceList
signal = cell(nr,1);
noise  = cell(nr,1);
idx    = cell(nr,1); 
bt     = cell(nr,1);

  
for ir = 1:nr
    
    % Peak amps
    [signal{ir},noise{ir},idx{ir}] = get_tdpa(trList,mRanges(ir,:),rRange,tRange,band);
         
    % Bin titles
    mbin      = floor(mRanges(ir,1));
    remainder = rem(mRanges(ir,1),1);
    if remainder==.5; annotTxt = 'high';
    else              annotTxt = 'low ';
    end
    bt{ir} = sprintf('M%i_{%s}',mbin,annotTxt);
end

ntr = cellfun(@(x) size(x,1), idx);     % No. of traces in each bin

% Find out whether it's broadband, strong motion or both
hasH = numel(find(strcmp('H',trList.instrCode)))>0;
hasL = numel(find(strcmp('L',trList.instrCode)))>0;
hasN = numel(find(strcmp('N',trList.instrCode)))>0;
hasG = numel(find(strcmp('G',trList.instrCode)))>0;
if (hasH  && ~hasL && ~hasN && ~hasG ); titleAppendix = 'BB'     ; end
if (~hasH && (hasL || hasN || hasG ) ); titleAppendix = 'SM'     ; end
if (hasH  && (hasL || hasN || hasG ) ); titleAppendix = 'SMandBB'; end



%% Estimate source durations
% ==========================

fprintf(1,'Computing source durations ... ')
m     = trList.m;
r     = trList.hypDist;
D_cdf = cell (nr,1);
t_cdf = cell (nr,1);
M_cdf = zeros(nr,1);


% b. Compute durations by summing up all sampled durations in each bin. This is equivalent to 
%    integrate out the magnitude.
for ir = 1:nr

    % Find median magnitudes in each bin
    idx_tmp = find( m>=mRanges(ir,1) & m<mRanges(ir,2) & r>=rRange(1) & r<rRange(2) );
    ntr_tmp = numel(idx_tmp);
    D_tmp   = zeros(ntr_tmp,nsample,'single');
    for i = 1:ntr_tmp
        D_tmp(i,:) = get_eqSourceDuration_WellsCoppersmith94(trList.m(idx_tmp(i)),nsample);
    end
    [D_cdf{ir},t_cdf{ir}] = ecdf(D_tmp(:));
end
fprintf(1,'done.\n')


%% KS-test & t-test for all mRange combinations
% =================

fprintf(1,'\nNOTE: KS-test is performed on log-amplitudes,\n      t-test  is performed on linear amplitudes\n')
% Compare distributions of consecutive pairs
tKS = zeros(nr,nr);
tST = zeros(nr,nr);

% Outer loop
for ir1 = 1:nr

    fprintf(1,sprintf('\n%i/%i\t',ir1,nr))

    % Linear and logarithmic amlitudes of two consecutive bins
    A1     = signal{ir1}.amps;
    A1_log = log10(A1);
    t1     = signal{ir1}.t;

    if ir1==4; 
        1+1;
    end
    
    for ir2 = 1: ir1
        A2     = signal{ir2}.amps;
        A2_log = log10(A2);
        t2     = signal{ir2}.t;

        nt = min([numel(t1), numel(t2)]);
        
        fprintf(1,sprintf('%i .. ',ir2))
        % (1) KS test can be done on linear amplitudes
        % (2) T-test assumes normal distribution, i.e. must be performed on
        %     log-amplitudes
        pks = zeros(nt,1);
        pst = zeros(nt,1);
        for it = 1:nt
            [~,pks(it)] = kstest2(A1    (:,it), A2    (:,it));
            [~,pst(it)] = ttest2 (A1_log(:,it), A2_log(:,it),'Vartype','unequal');
        end

        % Find time when pks and pst are last time above plevel
        1+1;
        %hold on; grid on; box on;
        %plot(t1,pks)
        %set(gca,'ylim',[-.1 1.1])
        idx_last=find(pks>=plevel,1,'last')+1;
        if isempty(idx_last); idx_last=1; end
        if idx_last>nt;      idx_last=nt; end
        %plot(t(idx_last),pks(idx_last),'or')
        tKS(ir1,ir2) = t1(idx_last);
    end
end

figure; clf; hold on; grid on; box on;
imagesc(tKS)
caxis([0 .1])
colorbar


1+1;
% %% Print figure
% %  ============
% if fig.print
%     set(gcf,'PaperPositionMode','auto')
%     %outPath = sprintf('~/programs/filterBank/fig/i%i/cf_tdpa/new/',iN);
%     
%     mString = strrep(sprintf('%3.1f_',fliplr(mRanges(:,1)')),'.','p');
%     %mString = sprintf('%i',fliplr(mRanges(:,1)'));
%     tString = sprintf('%i',tRange(2));
%     
%     if ~isnumeric(band)
%         %outFullName = sprintf('~/programs/filterBank/fig/i%i/cf_tdpa/new/%s_%s_%i_to_%ikm_%srec',iN,band,ornt,rRange,titleAppendix);
%         outFullName = sprintf('~/programs/filterBank/fig/i%i/cf_tdpa/new/%s_%s_%i_to_%ikm_m%s_t%ssec_%srec', ...
%             iN,band,ornt,rRange,mString,tString,titleAppendix);
%     else
%         outFullName = sprintf('~/programs/filterBank/fig/i%i/cf_tdpa/new/band_%i_%s_%i_to_%ikm_%srec',iN,band,ornt,rRange,titleAppendix);
%     end
%     
%     if strcmp(fig.xscale,'log'); outFullName = sprintf('%s_logt',outFullName); end
%     if strcmp(fig.xscale,'lin'); outFullName = sprintf('%s_lint',outFullName); end
%     
%     %print('-dpdf',[outFullName,'.pdf'])
%     print('-dpng','-r600',[outFullName,'.png'])
%     print('-depsc2',[outFullName,'.eps'])
% end