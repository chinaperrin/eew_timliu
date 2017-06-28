function summary_trList(TraceList,SkipList,opts)

global iN

% opts.rrange = [0 200];
% opts.dr     = 5;
% opts.mrange = [5 8];
% opts.dm     = .1;
% opts.mkSize = 4;
% summary_trList(zList,[],opts)

if nargin<3; 
    opts.plotData=0; 
    opts.printData=0; 
    opts.ftSize=12;
    opts.mkSize=3;
end

nmin = 100; % Min. no. of cases of skipped traces necessary
ftSize = opts.ftSize;
mkSize = opts.mkSize; 
mkSize_skipped = 6;

if opts.plotData
    figure(313); clf; 
    sp2=subplot(4,4,[4,8,12]);       hold on; grid on; box on;
    sp3=subplot(4,4,[13:15]);        hold on; grid on; box on;
    sp1=subplot(4,4,[1:3,5:7,9:11]); hold on; grid on; box on;
    %     L  = cell(50,1); 
    %     p  = cell(50,1);
    %ftSize = 15; 
end

fprintf(1,'\n--------------------------------------------------------------\n')

P      = [];
L      = {};
pcount = 0;

% TraceList ---------------------------------------------------------------
if ~isempty(TraceList)
    ntr = size(TraceList.fullName,1);
    fprintf(1,sprintf('%i traces were OK and made it into TraceList\n',ntr))
    
    nnocomment = sum(cellfun(@(x) isempty(x), TraceList.comment));
    fprintf(1,sprintf('  %i\t\t kept: no comment \n',nnocomment))
    
    % List commented traces
	commentList = {}; 
    for itr=1:ntr
        allComments = TraceList.comment{itr};
        commentList = unique([commentList; allComments]); 
    end
    
    ncmt   = numel(commentList);
    idxCmt = cell(ncmt,1);
    for icmt = 1:ncmt
        thisComment    = commentList{icmt};
        hasThisComment = cellfun(@(x) strcmp(thisComment,x), TraceList.comment,'uniformOutput',0);
        idxCmt{icmt}   = find(cellfun(@(x) sum(x), hasThisComment));
        ni             = numel(idxCmt{icmt});
        fprintf(1,sprintf('  %i\t\t kept: %s \n',ni,commentList{icmt}))
    end 

    
    if (opts.plotData)
    
        idx_scsn = find( strcmp(TraceList.dataSetName,'scsn')      |strcmp(TraceList.dataSetName,'scsnPx') );
        idx_jp   = find( strcmp(TraceList.dataSetName,'kNet')      |strcmp(TraceList.dataSetName,'kikNet') |strcmp(TraceList.dataSetName,'kNet_manPx') |strcmp(TraceList.dataSetName,'kikNet_manPx') ); 
        idx_nga  = find( strcmp(TraceList.dataSetName ,'ngawest1') |strcmp(TraceList.dataSetName ,'ngawest1_manPx'));
        idx_wch  = find( strcmp(TraceList.dataSetName ,'wenchuan'));

        if ~isempty(idx_jp);
            p  = plot(TraceList.dist.hyp(idx_jp)  ,TraceList.eq.m(idx_jp)  ,'or','markerFaceColor','r','markerSize',mkSize);
            pcount = pcount+1;
            P      = [P; p];
            L      = [L,'K-NET & KiK-net'];
        end
        
        if ~isempty(idx_scsn);
            p      = plot(TraceList.dist.hyp(idx_scsn),TraceList.eq.m(idx_scsn),'^b','markerFaceColor','b','markerSize',mkSize);
            pcount = pcount+1;
            P      = [P; p];
            L      = [L,'SCSN'];
        end
        
        if ~isempty(idx_nga);  
            p = plot(TraceList.dist.hyp(idx_nga) ,TraceList.eq.m(idx_nga) ,'sk','markerFaceColor','k','markerSize',mkSize); 
        	pcount = pcount+1;
            P      = [P; p];
            L      = [L,'NGA West 1'];
        end
        
        if ~isempty(idx_wch);
            p = plot(TraceList.dist.hyp(idx_wch) ,TraceList.eq.m(idx_wch) ,'dk','markerFaceColor','c','markerSize',mkSize); 
        	pcount = pcount+1;
            P      = [P; p];
            L      = [L,'Wenchuan 2008'];
        end
    end
end


% SkipList ----------------------------------------------------------------
if ~isempty(SkipList)
    nsk = size(SkipList.fullName,1);
    fprintf(1,sprintf('\n%i traces were SKIPPED\n',nsk))

    %     % Here's an alternative of how you can scan for a particular comment
    %     tmp           = cellfun(@(x) regexp(x,'has unclipped BB record'), TraceList.comment,'uniformOutput',0);
    %     idxHasComment = find(cellfun(@(x) sum(cell2mat(x)), tmp));
    %     idxNotComment  = find(cellfun(@(x) sum(cell2mat(x)), tmp)==0);

    skipReasonList = {}; 
    for isk=1:nsk;
        allComments    = SkipList.comment{isk};
        skipReasonList = unique([skipReasonList; allComments]); 
    end
    nsr = numel(skipReasonList);
    
    idx = cell(nsr,1);
    for isr = 1:nsr 
        thisReason    = skipReasonList{isr};
        hasThisReason = cellfun(@(x) strcmp(thisReason,x), SkipList.comment,'uniformOutput',0);
        idx{isr}      = find(cellfun(@(x) sum(x), hasThisReason));
        ni            = numel(idx{isr});
        fprintf(1,sprintf('  %i\t\t skipped: %s \n',ni,skipReasonList{isr}))
    end 
    
    if (opts.plotData)
        
        symbolList = {'*','v','^','p','d','s','o','h','<','>','.'};
        ns         = numel(symbolList);
        pcount  = 0;
        for isr = 1:nsr
            
            ni = numel(idx{isr});
            %fprintf(1,sprintf('  %i\t\t skipped: %s \n',ni,skipReasonList{isr}))
    
            if ni>nmin
                pcount = pcount+1;
                if pcount<=ns; symbol = symbolList{pcount};
                               sColor = [0 .4 0];
                else           symbol = symbolList{pcount-ns};
                               sColor = 'm';
                end
                p = plot(SkipList.dist.hyp(idx{isr}), SkipList.eq.m(idx{isr}), symbol,'markerSize',mkSize_skipped,'color',sColor);
                P      = [P; p];
                L      = [L,skipReasonList{isr}];
            end
        end
    end
end
fprintf(1,'\n  Note: each trace can have multiple comments, \n  some of which warrant skipping, some of which do not.\n')


if opts.plotData     
    % Legend
    %     firstIdx = find(cellfun(@(x) ~isempty(x), p),1,'first');
    %     lastIdx  = find(cellfun(@(x) ~isempty(x), p),1,'last');
    %     L        = L(firstIdx:lastIdx);
    %     p        = p(firstIdx:lastIdx);
    %l1       = legend(cell2mat(P),L);
    l1       = legend(P,L);
    set(l1,'fontSize',ftSize-2,'location','northWest')
    
    %set(gca,'xlim',[-1.5 50],'ylim',[2.9 9.1])
    %set(gca,'xlim',[-1.5 100],'ylim',[1.9 8.1],'fontSize',ftSize)
    set(gca,'xlim',[-1.5 ceil(max(TraceList.dist.hyp))+1], ...
        'ylim', [floor(min(TraceList.eq.m))-0.1 ceil(max(TraceList.eq.m))],'fontSize',ftSize)
    xlabel('Hypocentral Distance [km]','fontSize',ftSize)
    ylabel('Catalog Magnitude','fontSize',ftSize)
    set(gca,'xlim',opts.rrange,'ylim',opts.mrange,'fontSize',ftSize)
    
    % Title
    tString = sprintf('No. of records: kNet & kikNet (%i), NGA West 1 (%i), SCSN (%i)',numel(idx_jp),numel(idx_nga),numel(idx_scsn));
    fprintf(1,'not adding title\n')
    %title(tString,'fontSize',ftSize)
    
    
    subplot(sp3);
    hist(TraceList.dist.hyp,opts.rrange(1):opts.dr:opts.rrange(2));
    xlabel('Hypocentral Distance [km]','fontSize',ftSize) 
    ylabel('No. of cases','fontSize',ftSize)
    set(gca,'xlim',opts.rrange,'fontSize',ftSize)
    set(get(gca,'child'),'FaceColor','k','EdgeColor','w');
    
    subplot(sp2);
    hist(TraceList.eq.m,opts.mrange(1):opts.dm:opts.mrange(2));
    xlabel('Catalog Magnitude','fontSize',ftSize)
    ylabel('No. of cases','fontSize',ftSize)
    set(gca,'xlim',opts.mrange,'yscale','lin','fontSize',ftSize)
    set(get(gca,'child'),'FaceColor','k','EdgeColor','w');
    camroll(90)
    

    if opts.printData
        figFullName = '~/programs/seismo/fig/i39/data/new/data';
        set(gcf,'PaperPositionMode','auto')
        %print('-depsc2',[figName,'.eps'])
        print('-dpng',figFullName)
    end
end