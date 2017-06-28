function summary_trList(TraceList,SkipList,o_plotData,o_printData)

global iN

if nargin<3; o_plotData  = 0; end
if nargin<4; o_printData = 0; end

nmin = 100; % Min. no. of cases of skipped traces necessary
mkSize = 3;
mkSize_skipped = 6;

if (o_plotData)
    %figure(313); clf; hold on; grid on; box on;
%     L  = cell(50,1); 
%     p  = cell(50,1);
    
    ftSize = 15; 
end

fprintf(1,'\n--------------------------------------------------------------\n')

P      = [];
L      = {};
pcount = 0;   

% TraceList ---------------------------------------------------------------
if ~isempty(TraceList)
    ntr = size(TraceList.fullName,1);
    fprintf(1,sprintf('%i\ttraces were ok and made it into TraceList\n',ntr))
    
    if (o_plotData)
    
        idx_scsn = find( strcmp(TraceList.dataSetName,'scsn')      |strcmp(TraceList.dataSetName,'scsnPx') );
        idx_jp   = find( strcmp(TraceList.dataSetName,'kNet')      |strcmp(TraceList.dataSetName,'kikNet') |strcmp(TraceList.dataSetName,'kNet_manPx') |strcmp(TraceList.dataSetName,'kikNet_manPx') ); 
        idx_nga  = find( strcmp(TraceList.dataSetName ,'ngawest1') |strcmp(TraceList.dataSetName ,'ngawest1_manPx'));
        idx_wch  = find( strcmp(TraceList.dataSetName ,'wenchuan'));

        if ~isempty(idx_scsn); 
            p      = plot(TraceList.hypDist(idx_scsn),TraceList.m(idx_scsn),'^b','markerFaceColor','b','markerSize',mkSize);
            pcount = pcount+1;
            P      = [P; p];
            L      = [L,'SCSN'];
        end
        
        if ~isempty(idx_jp);
            p  = plot(TraceList.hypDist(idx_jp)  ,TraceList.m(idx_jp)  ,'or','markerFaceColor','r','markerSize',mkSize);
            pcount = pcount+1;
            P      = [P; p];
            L      = [L,'kNet & kikNet'];
        end
        
        if ~isempty(idx_nga);  
            p = plot(TraceList.hypDist(idx_nga) ,TraceList.m(idx_nga) ,'sk','markerFaceColor','k','markerSize',mkSize); 
        	pcount = pcount+1;
            P      = [P; p];
            L      = [L,'NGA West 1'];
        end
        
        if ~isempty(idx_wch);
            p = plot(TraceList.hypDist(idx_wch) ,TraceList.m(idx_wch) ,'dk','markerFaceColor','c','markerSize',mkSize); 
        	pcount = pcount+1;
            P      = [P; p];
            L      = [L,'Wenchuan 2008'];
        end
    end
end


% SkipList ----------------------------------------------------------------
if ~isempty(SkipList)
    nsk = size(SkipList.fullName,1);
    fprintf(1,sprintf('%i\ttraces were skipped\n\n',nsk))
    
    tmpList = {}; 
    for isk=1:nsk;
        allComments = SkipList.comment(isk);
        if numel(allComments)>1;
            1+1;
        end
        tmpList = [tmpList; allComments]; 
    end
    skpReasonList = unique(tmpList);
    nsr           = numel(skpReasonList);
    %hasComment=cell2mat(cellfun(@(x) ~isempty(x) , SkipList.comment,'uniformOutput',0));
    %skpReasonList = unique(SkipList.comment(hasComment));
    
    idx = cell(nsr,1);
    for isr = 1:nsr 
        thisReason    = skpReasonList{isr};
        hasThisReason = cellfun(@(x) strcmp(thisReason,x), SkipList.comment,'uniformOutput',0);
        idx{isr}      = find(cellfun(@(x) sum(x), hasThisReason));
        ni            = numel(idx{isr});
        %hasThisReason = cellfun(@(x) ~isempty(regexp(thisReason,x)), SkipList.comment,'uniformOutput',0);
        %idx{isr}      = find(cell2mat(cellfun(@(x) ~isempty(x{1}), hasThisReason,'uniformOutput',0)));
        %idx{isr}      = find(cell2mat(cellfun(@(x) sum(x), hasThisReason,'uniformOutput',0)));
        %idx{isr} = find(strcmp(SkipList.comment,reason));   
        fprintf(1,sprintf('  %i\t\t skipped: %s \n',ni,skpReasonList{isr}))
    end 
    
    if (o_plotData)
        
        symbolList = {'*','v','^','p','d','s','o','h','<','>','.'};
        ns         = numel(symbolList);
        pcount  = 0;
        for isr = 1:nsr
            
            ni = numel(idx{isr});
            %fprintf(1,sprintf('  %i\t\t skipped: %s \n',ni,skpReasonList{isr}))
    
            if ni>nmin
                pcount = pcount+1;
                if pcount<=ns; symbol = symbolList{pcount};
                               sColor = [0 .4 0];
                else           symbol = symbolList{pcount-ns};
                               sColor = 'm';
                end
                p = plot(SkipList.hypDist(idx{isr}), SkipList.m(idx{isr}), symbol,'markerSize',mkSize_skipped,'color',sColor);
                P      = [P; p];
                L      = [L,skpReasonList{isr}];
            end
        end
    end
end

if o_plotData 
    
    % Legend
    %     firstIdx = find(cellfun(@(x) ~isempty(x), p),1,'first');
    %     lastIdx  = find(cellfun(@(x) ~isempty(x), p),1,'last');
    %     L        = L(firstIdx:lastIdx);
    %     p        = p(firstIdx:lastIdx);
    %l1       = legend(cell2mat(P),L);
    l1       = legend(P,L);
    set(l1,'fontSize',ftSize,'location','northWest')
    
    %set(gca,'xlim',[-1.5 50],'ylim',[2.9 9.1])
    %set(gca,'xlim',[-1.5 100],'ylim',[1.9 8.1],'fontSize',ftSize)
    set(gca,'xlim',[-1.5 ceil(max(TraceList.hypDist))+1], ...
        'ylim', [floor(min(TraceList.m))-0.1 ceil(max(TraceList.m))],'fontSize',ftSize)
    xlabel('Hypocentral Distance [km]','fontSize',ftSize)
    ylabel('Catalog Magnitude','fontSize',ftSize)
    set(gca,'fontSize',ftSize)
    
    % Title
    tString = sprintf('No. of records: kNet & kikNet (%i), NGA West 1 (%i), SCSN (%i)',numel(idx_jp),numel(idx_nga),numel(idx_scsn));
    fprintf(1,'not adding title\n')
    %title(tString,'fontSize',ftSize)
    
    if o_printData
        figDir  = strcat(['~/programs/filterBank/fig/i',num2str(iN),'/data/new/']);
        figName = sprintf('%ssummary_%itraces',figDir,ntr);
        set(gcf,'PaperPositionMode','auto')
        %set(gcf,'PaperPositionMode','manual')
        %set(gcf,'DefaultFigureColor','remove')
        %set(gcf,'DefaultFigureColor',[1 1 1])
        
        print('-depsc2',[figName,'.eps'])
        print('-dpng',[figName,'.png'])
    end
end