function summary_trList(TraceList,SkipList,o_plotData,o_printData)

global iN

if nargin<3; o_plotData  = 0; end
if nargin<4; o_printData = 0; end

if (o_plotData)
    figure(313); clf; hold on; grid on; box on;
    L  = cell (50,1); 
    p  = zeros(50,1);
    
    ftSize = 15; 
end

fprintf(1,'\n--------------------------------------------------------------\n')

% TraceList ---------------------------------------------------------------
if ~isempty(TraceList)
    ntr = size(TraceList.fullName,1);
    fprintf(1,sprintf('%i\ttraces were ok and made it into TraceList\n',ntr))
    
    if (o_plotData)
    
        idx_scsn = find( strcmp(TraceList.dataSetName,'scsn') |strcmp(TraceList.dataSetName,'scsnPx') );
        idx_jp   = find( strcmp(TraceList.dataSetName,'kNet') |strcmp(TraceList.dataSetName,'kikNet') );
        idx_nga  = find(strcmp(TraceList.dataSetName,'ngawest1'));

        p(3) = plot(TraceList.hypDist(idx_scsn),TraceList.m(idx_scsn),'^b','markerFaceColor','b','markerSize',5);
        p(1) = plot(TraceList.hypDist(idx_jp)  ,TraceList.m(idx_jp)  ,'or','markerFaceColor','r','markerSize',5);
        p(2) = plot(TraceList.hypDist(idx_nga) ,TraceList.m(idx_nga) ,'sk','markerFaceColor','k','markerSize',5);
        L{3} = 'SCSN';     
        L{1} = 'kNet & kikNet';
        L{2} = 'NGA West 1';
        
    end
end


% SkipList ----------------------------------------------------------------
if ~isempty(SkipList)
    nsk = size(SkipList.fullName,1);
    fprintf(1,sprintf('%i\ttraces were skipped\n\n',nsk))
    
    skpReasonList = unique(SkipList.comment);
    nsr           = numel(skpReasonList);
    
    idx = cell(nsr,1);
    
    for i = 1:nsr 
        reason = skpReasonList{i};
        idx{i} = find(strcmp(SkipList.comment,reason));
            
        ni = numel(idx{i});
        fprintf(1,sprintf('  %i\t\t skipped: %s \n',ni,skpReasonList{i}))
    end
    
    if (o_plotData)
        
        symbolList = {'*','v','^','p','^','d','s','o','h','<','>','.'};
        plotCount  = 0;
        for i = 1:nsr
            
            ni = numel(idx{i});
            %fprintf(1,sprintf('  %i\t\t skipped: %s \n',ni,skpReasonList{i}))
    
            if ni>2000
                plotCount      = plotCount+1;
                symbol         = symbolList{plotCount};
                p(plotCount+3) = plot(SkipList.hypDist(idx{i}), SkipList.m(idx{i}), symbol,'markerSize',3);
                L{plotCount+3} = skpReasonList{i};
            end
        end
    end
end

if o_plotData
    
    % Legend
    firstIdx = find(cellfun(@(x) ~isempty(x), L),1,'first');
    lastIdx  = find(cellfun(@(x) ~isempty(x), L),1,'last');
    L        = L(firstIdx:lastIdx);
    p        = p(firstIdx:lastIdx);
    l1       = legend(p,L);
    set(l1,'fontSize',ftSize,'location','southEast')
    
    %set(gca,'xlim',[-1.5 50],'ylim',[2.9 9.1])
    set(gca,'xlim',[-1.5 100],'ylim',[1.9 8.1],'fontSize',ftSize)
    xlabel('Hypocentral Distance [km]','fontSize',ftSize)
    ylabel('Catalog Magnitude','fontSize',ftSize)
    set(gca,'fontSize',ftSize)
    
    % Title
    tString = sprintf('No. of records: kNet & kikNet (%i), NGA West 1 (%i), SCSN (%i)',numel(idx_jp),numel(idx_nga),numel(idx_scsn));
    title(tString,'fontSize',ftSize)
    
    if o_printData
        figDir  = strcat(['~/programs/filterBank/fig/i33/data/new/']);
        figName = sprintf('%ssummary_%itraces',figDir,ntr);
        set(gcf,'PaperPositionMode','auto')
        %set(gcf,'PaperPositionMode','manual')
        %set(gcf,'DefaultFigureColor','remove')
        %set(gcf,'DefaultFigureColor',[1 1 1])
        
        print('-depsc2',[figName,'.eps'])
        print('-dpng',[figName,'.png'])


    end
end