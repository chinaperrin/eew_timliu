% Plot attenuation curves for different magnitude bins, in order to find
% erroneous data.
%
% menandrin@gmail.com, 130624


%% 
addpath(genpath('../../../../matlab/fct_downloads/'))
addpath(genpath('../../../../../gmpe/'))
addpath(genpath('../../../../../matlab/sac/'))
addpath(genpath('../../../../../matlab/myFct/maps/'))
addpath(genpath('../../../../../VS/data/plot/'))
clear all



%% ===============
%  0. Preparations

mmin = 5;
mmax = 7.5;
dm   = 0.5;
mm   = mmin:dm:mmax;
nm   = numel(mm)-1;

xband = 4;
xsnp  = 6;

alogmax = 0;
alogmin = -6;

level = 2;         % Level for percentiles

% Prepare distribution arrays
dint   = 20;
dedges = 0:dint:220;
nedges = numel(dedges)-1;
dd     = 0:220;



%%
% Station list
fprintf(1,'\n\n\n--------------------------------------------------------------------------------------')
fprintf(1,'\nImporting station lists\n')

% SoCal
stFileName         = '../../../../data/stations/eewvs-stations_allCali.hinv.fmtd2';
[stationList_cali] = import_stationlist_fB(stFileName);

stFileName     = '../../../../data/stations/japan/sitepub_all_en.txt';
stationList_jp = import_k_kik_net_list(stFileName);




%% ==================
%  1.1 Load data sets
%  ==================
fprintf(1,'\n\nLoading filterBank.m-output ...   \n')

% All
% dataSetNames = {'/scratch/memeier/VS/data/wform/socal/scsn_900101_011231/9bands/'; ...
%                 '/scratch/memeier/VS/data/wform/socal/scsn_020101_031231/9bands/'; ...
%                 '/scratch/memeier/VS/data/wform/socal/scsn_040101_051231/9bands/'; ...
%                 '/scratch/memeier/VS/data/wform/socal/scsn_060101_071231/9bands/'; ...
%                 '/scratch/memeier/VS/data/wform/socal/scsn_080101_091231/9bands/'; ...
%                 '/scratch/memeier/VS/data/wform/socal/scsn_100101_111231/9bands/'; ...
%                 '/scratch/memeier/VS/data/wform/socal/scsn_120101_131231/9bands/'};

dataSetNames = {'/scratch/memeier/VS/data/wform/socal/scsn_900101_011231/9bands_i24/'};

[traceList,skipList,fc,snpLength,nsnp] = import_fbOut(dataSetNames);

% nnn=size(traceList.m,1)-size(traceList.snr,1);
% traceList.snr = [traceList.snr; zeros(nnn,1)];
% nnn=size(skipList.m,1)-size(skipList.snr,1);
% skipList.snr = [skipList.snr; zeros(nnn,1)];

% Extract only horizontal components
tIdx  = find(logical(strcmp(traceList.orntCode,'E')) | logical(strcmp(traceList.orntCode,'N')));
tList = traceList.selectSubList(tIdx);
sIdx  = find(logical(strcmp(skipList.orntCode,'E')) | logical(strcmp(skipList.orntCode,'N')));
sList = skipList.selectSubList(sIdx);


%%
tL_idx_up = [];
tL_idx_lo = [];
drel_outUp   = []; 
drel_outLo   = []; 

for im = 1:nm
    
    m_lo = mm(im);
    m_up = mm(im+1);
    idx  = find( (tList.m>=m_lo) & (tList.m<m_up) );
    idx2 = find( (sList.m>=m_lo) & (sList.m<m_up) );
    
    % Prepare Plot
    figure(342); clf; hold on;
    grid on; %set(gca,'xscale','log')
    xlabel('Hypocentral distance [km]','fontSize',15)
    ylabel('log_{10}(PGV) [m/s]','fontSize',15)
    title(['Magnitude range ', ...
        num2str(m_lo),' <= M < ',num2str(m_up)],'fontSize',15)
    set(gca,'xlim',[0.1 220],'fontSize',15)
    set(gca,'ylim',[alogmin alogmax],'fontSize',15)
    set(gca,'ytick',alogmin:1:alogmax)
    whitebg(gcf,'k')
    set(gcf,'InvertHardcopy','off')
      
    
    % Socal   -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    if (~isempty(idx))

        dVal = tList.hypDist(idx);
        aVal = tList.pgv(idx);
        ornt = tList.orntCode{idx};
        h1   = plot(dVal,log10(aVal),'xw','lineWidth',2);
        
%         % Locate, plot and save outliers
%         [~,up,~,~] = CH2007_2sigma(m_up,dVal,'PGV','H','S','S',0);
%         up         = up/1e2;
%         diff       = aVal-up;
%         idx_up     = find(diff>0);              
%         tL_idx_up  = [tL_idx_up; idx(idx_up)];
%         dtmp       = aVal(idx_up)./up(idx_up);    
%         drel_outUp = [drel_outUp; dtmp];
% 
%         [~,~,lo,~] = CH2007_2sigma(m_lo,dVal,'PGV','H','S','R',0);
%         lo         = lo/1e2;
%         diff       = -aVal+lo;
%         idx_lo     = find(diff>0);              
%         tL_idx_lo  = [tL_idx_lo; idx(idx_lo)];
%         dtmp       = aVal(idx_lo)./lo(idx_lo);    
%         drel_outLo = [drel_outLo; dtmp];
        
        
    end
     
         % Socal   -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    if (~isempty(idx2))

        dVal2 = sList.hypDist(idx2);
        aVal2 = sList.pgv(idx2);
        h2    = plot(dVal2,log10(aVal2),'xy','lineWidth',2);
        
    end
     
     
     % Compute, transform and plot GMPEs
     [Y,~,~,~]   = CH2007((m_up+m_lo)/2,dd,'PGV','H','S','R',0);
     [~,UP,~,~]  = CH2007( m_up        ,dd,'PGV','H','S','S',0);
     [~,~,LO,~]  = CH2007( m_lo        ,dd,'PGV','H','S','R',0);
     [~,UP2,~,~] = CH2007_2sigma( m_up ,dd,'PGV','H','S','S',0);
     [~,~,LO2,~] = CH2007_2sigma( m_lo ,dd,'PGV','H','S','R',0);
     
     Y   = Y/1e2;                                       % Convert to m/s
     UP  = UP/1e2;
     UP2 = UP2/1e2;
     LO  = LO/1e2;
     LO2 = LO2/1e2;
     
     h4  = plot(dd,log10(Y)  ,'-r','lineWidth',2);
     h4b = plot(dd,log10(UP) ,'--r');
     plot(dd,log10(UP2),'--r');
     h4c = plot(dd,log10(LO) ,':r','lineWidth',2);
     plot(dd,log10(LO2),':r','lineWidth',2);
     legend([h1,h2,h4],'accepted','''outlier''','CH2007 (M5.25)')
     %print('-dpng','-r150',['gmpe_' num2str(im) '.png']);
     
     1+1;
     pause
     
end










%% APPENDIX
% %% Add outliers to separate traceList and sort them
% 
% % SoCal - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% if (exist('idx','var'))
%     
%     % Create sublist with outliers and sort such that largest errors come first
%     OutUp                    = tList.selectSubList(tL_idx_up);
%     [drel_outUp_std,idx_std] = sort(drel_outUp,'descend');
%     OutUp.sortList(idx_std);
%     
%     OutLo                    = tList.selectSubList(tL_idx_lo);
%     [drel_outLo_std,idx_std] = sort(drel_outLo,'ascend');
%     OutLo.sortList(idx_std);
%     
%     % Plot errors   
%     figure(33);clf; hold on
%     e1 = plot(drel_outLo_std,'ob');
%     e2 = plot(drel_outUp_std,'xk');
%     set(gca,'yscale','log')
%     grid on
%     l1 = legend([e1; e2],'Outliers << CH2007','Outliers >> CH2007','fonzSize',15);
%     set(l1,'fontSize',15)
%     ylabel('relative amplitude wrt CH2007 (aVal/CH07)','fontSize',15)
% end
% 
% 
% 
% % SoCal90ies  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% if (exist('idx90','var'))
% 
%     % Create sublist with outliers and sort such that largest errors come first
%     OutUp_sc90                      = tList90.selectSubList(tL_idx_up_sc90);
%     [drel_outUp_sc90_std,sortIndex] = sort(drel_outUp_sc90,'descend');
%     OutUp_sc90.sortList(sortIndex);
%     
%     OutLo_sc90                      = tList90.selectSubList(tL_idx_lo_sc90);
%     [drel_outLo_sc90_std,sortIndex] = sort(drel_outLo_sc90,'ascend');
%     OutLo_sc90.sortList(sortIndex);
%     
%     % Plot errors   - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%     figure(34);clf; hold on
%     e1 = plot(drel_outLo_sc90_std,'ob');
%     e2 = plot(drel_outUp_sc90_std,'xk');
%     set(gca,'yscale','log')
%     grid on
%     l1 = legend([e1; e2],'Outliers << CH2007','Outliers >> CH2007','fonzSize',15);
%     set(l1,'fontSize',15)
%     ylabel('relative amplitude wrt CH2007 (aVal/CH07)','fontSize',15)
% 
% end
% 
% 
% 
% % Japan       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 
% 
% 
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%
% %  Analyse oulier lists  %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % Preselection, e.g. to get rid of traces with only small errors
% %OutLo_mod = OutLo.selectSubList(1:4000);
% 
% 
% % Choose traceLists. Analyse all entries of <smallList> and use <largeList>
% % to find e.g. co-recorded traces or traces from the same station.
% largeList = tList;
% %largeList = tList90;
% 
% %smallList = tList90;
% %smallList = OutLo_sc90;
% %smallList = OutUp_sc90;
% smallList = OutUp;
% %smallList = OutLo_mod;
% 
% 
% % A. Analyse station lists of outliers
% % ====================================
% n       = numel(smallList.sRate);
% staList = cell(n,1);
% 
% for i = 1:n
%     [~,~,~,sta_tmp,~] = get_stationCoords_cali(smallList.fullName{i},stationList_cali);
%     staList{i}        = sta_tmp;
% end
% [countList] = count_station_instances(staList);
% 
% nsta = size(countList,1);
% for ista = 1:nsta
%    
%     fprintf(1,['Looking at ',num2str(ista),'th of ',num2str(nsta),' stations \n'])
%    
%    stationName = countList{ista,1};
%    [idx]        = find_station_traces(stationName,largeList);
%    
%    nout              = countList{ista,2};
%    ntot              = numel(idx);
%    countList{ista,3} = nout/ntot*100;       % Add fraction of outliers to countList
% end
% 
% 
% 
% % B. Go through individual outlier traces
% % =======================================
% 
% ntr = numel(smallList.sRate);
% for itr = 1:ntr
%     
%     traceFullName = smallList.fullName{itr};
%     slashIdx      = regexp(traceFullName,'/');
%     recordName    = traceFullName(slashIdx(end)+1:end);
%     fprintf(1,['\n\n\n     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n'])
%     fprintf(1,['     Looking at ',recordName,'\n'])
%     fprintf(1,['     . . . . . . . . . . . . . . . . .\n\n'])
%     
%     
%     
%     % Plot all amplitudes from the same station
%     % -----------------------------------------
%     fprintf(1,['Searching for traces from the same station ...\n'])
%     
%     % Make traceList with all traces of a station
%     stationName              = recordName(end-10:end-8);
%     [traceIndices]           = find_station_traces(stationName,largeList);
%     trL_sameStation          = largeList.selectSubList(traceIndices);
%     
%     [ort,~] = find_orientation(trL_sameStation);
%     lgc_z = cellfun(@(x) strcmp(x,'Z'), ort);
%     lgc_n = cellfun(@(x) strcmp(x,'N'), ort);
%     lgc_e = cellfun(@(x) strcmp(x,'E'), ort);
%     
%     figure(12); clf; hold on
%     d = trL_sameStation.hypDist;
%     a = log10(cell2mat(trL_sameStation.amax));
%     m = trL_sameStation.m;
%     scatter(d(lgc_z),a(lgc_z),50,m(lgc_z),'filled','marker','d')
%     scatter(d(lgc_n),a(lgc_n),120,m(lgc_n),'filled','marker','p')
%     scatter(d(lgc_e),a(lgc_e),50,m(lgc_e),'filled','marker','^')
%     
%     pd = plot(smallList.hypDist(itr),log10(smallList.amax{itr}), ...
%         'pk','markerSize',14,'markerFaceColor','r');
%     
%     grid on
%     colorbar
%     title(['All traces of station ',stationName],'fontSize',15)
%     set(gca,'xscale','log','xlim',[0.1 1000],'ylim',[alogmin alogmax],'ytick', ...
%         alogmin:1:alogmax,'clim',[3 8])
%     
%     % GMPEs
%     magnitude  = smallList.m(itr);
%     [Ym,~,~,~] = CH2007(magnitude,dd,'PGV','H','S','R',0);
%     pm         = plot(dd,log10(Ym/1e2),'-r','lineWidth',2);
%     
%     [Y3,~,~,~]  = CH2007(3,dd,'PGV','H','S','R',0);
%     [Y4,~,~,~]  = CH2007(4,dd,'PGV','H','S','R',0);
%     [Y5,~,~,~]  = CH2007(5,dd,'PGV','H','S','R',0);
%     [Y6,~,~,~]  = CH2007(6,dd,'PGV','H','S','R',0);
%     [Y7,~,~,~]  = CH2007(7,dd,'PGV','H','S','R',0);
%     p3 = plot(dd,log10(Y3/1e2),'-k');
%     p4 = plot(dd,log10(Y4/1e2),'--k');
%     p5 = plot(dd,log10(Y5/1e2),':k');
%     p6 = plot(dd,log10(Y6/1e2),'-b');
%     p7 = plot(dd,log10(Y7/1e2),'--b');
%     l1 = legend([pd;pm;p7;p6;p5;p4;p3],'Current data point',['CH2007 (M',num2str(magnitude),')'],'M7','M6','M5','M4','M3');
%     set(l1,'location','SouthWest')
%     
%     
%     
%     % Plot all amplitudes of same event
%     % ---------------------------------
%     fprintf(1,['Searching for traces from the same event ...\n'])
% 
%     ptIdx         = regexp(recordName,'\.');
%     eventName     = recordName(1:ptIdx(1)-1);
% 
%     [traceIndices] = find_event_traces(eventName,largeList);
%     trL_sameEvent  = largeList.selectSubList(traceIndices);
%     
%     figure(13); clf; hold on
%     d = trL_sameEvent.hypDist;
%     a = log10(cell2mat(trL_sameEvent.amax));
%     m = trL_sameEvent.m;
%     scatter(d,a,30,m,'filled')
%     pd = plot(smallList.hypDist(itr),log10(smallList.amax{itr}), ...
%         'pk','markerSize',14,'markerFaceColor','r');
%     
%     grid on
%     colorbar
%     title(['All traces of event ',eventName],'fontSize',15)
%     set(gca,'xscale','log','xlim',[0.1 1000],'ylim',[alogmin alogmax],'ytick', ...
%         alogmin:1:alogmax,'clim',[3 8])
%     
%     % GMPEs
%     magnitude  = smallList.m(itr);
%     [Ym,~,~,~] = CH2007(magnitude,dd,'PGV','H','S','R',0);
%     pm         = plot(dd,log10(Ym/1e2),'-r','lineWidth',2);
%     
%     [Y3,~,~,~]  = CH2007(3,dd,'PGV','H','S','R',0);
%     [Y4,~,~,~]  = CH2007(4,dd,'PGV','H','S','R',0);
%     [Y5,~,~,~]  = CH2007(5,dd,'PGV','H','S','R',0);
%     [Y6,~,~,~]  = CH2007(6,dd,'PGV','H','S','R',0);
%     [Y7,~,~,~]  = CH2007(7,dd,'PGV','H','S','R',0);
%     p3 = plot(dd,log10(Y3/1e2),'-k');
%     p4 = plot(dd,log10(Y4/1e2),'--k');
%     p5 = plot(dd,log10(Y5/1e2),':k');
%     p6 = plot(dd,log10(Y6/1e2),'-b');
%     p7 = plot(dd,log10(Y7/1e2),'--b');
%     l2 = legend([pd;pm;p7;p6;p5;p4;p3],'Current data point',['CH2007 (M',num2str(magnitude),')'],'M7','M6','M5','M4','M3');
%     set(l2,'location','SouthWest')
% 
%     
%     % Plot waveforms
%     % ------------------------
%     plotTraceSummary(traceFullName,largeList,60,1,10,1)
% 
% %    if (strcmp(recordName(end-4),'Z'))
%         1+1;
% %    end
% end
% 
% 
% 
% 
% 
% %%
% 
% 
% 
% %n       = 3000;
% n       = numel(OutLo.sRate);
% staList = cell(n,1);
% 
% for i = 1:n
%     %[~,~,~,sta_tmp,~] = get_stationCoords_cali(OutUp.fullName{i},stationList_cali);
%     [~,~,~,sta_tmp,~] = get_stationCoords_cali(OutLo.fullName{i},stationList_cali);
%     staList{i}        = sta_tmp;
% end
% [countList] = count_station_instances(staList);
% 
% 
% 
% nsta = 10;
% for i = 1:nsta
%     
%     % Make traceList with all traces of a station
%     [traceIndices]           = find_station_traces(countList{i},tList);
%     flg_select               = false(numel(tList.t0),1);
%     flg_select(traceIndices) = true;
%     trL_tmp                  = tList.selectSubList(flg_select);
%     
%     figure(12); clf; hold on
%     d = trL_tmp.hypDist;
%     a = log10(cell2mat(trL_tmp.amax));
%     m = trL_tmp.m;
%     scatter(d,a,10,m,'filled')
%     set(gca,'xscale','log')
%     set(gca,'xlim',[0.1 1000])
%     %set(gca,'ylim',[alogmin alogmax])
%     set(gca,'ytick',[alogmin:1:alogmax])
%     grid on
%     colorbar
%     title(['Station ',countList{i,1},' (',num2str(countList{i,2}),' outliers)'],'fontSize',15)
%     
%     % GMPEs
%     [Y3,~,~,~]  = CH2007(3,dd,'PGV','H','S','R',0);
%     [Y4,~,~,~]  = CH2007(4,dd,'PGV','H','S','R',0);
%     [Y5,~,~,~]  = CH2007(5,dd,'PGV','H','S','R',0);
%     [Y6,~,~,~]  = CH2007(6,dd,'PGV','H','S','R',0);
%     [Y7,~,~,~]  = CH2007(7,dd,'PGV','H','S','R',0);
%     p3=plot(dd,log10(Y3/1e2),'-k');
%     p4=plot(dd,log10(Y4/1e2),'--k');
%     p5=plot(dd,log10(Y5/1e2),':k');
%     p6=plot(dd,log10(Y6/1e2),'-b');
%     p7=plot(dd,log10(Y7/1e2),'--b');
%     legend([p3;p4;p5;p6;p7],'M3','M4','M5','M6','M7')
%     pause
% end
% 
% 
% idx=find(cell2mat(trL_tmp.amax)<5e-7);
% plot(trL_tmp.hypDist(idx),log10(cell2mat(trL_tmp.amax(idx))),'xr');
% 
% 
% % Traces with amax = 0  - - - - - - - - - - - - - - - - - - - - - - - - - -  
% idx_zilt   = find(cell2mat(OutLo.amax)==0);
% OutLo_zilt = OutLo.selectSubList(idx_zilt);
% 
% n       = numel(OutLo_zilt.t0);
% staList = cell(n,1);
% for i = 1:n
%     [~,~,~,sta_tmp,~] = get_stationCoords_cali(OutLo_zilt.fullName{i},stationList_cali);
%     staList{i}        = sta_tmp;
% end
% [countList] = count_station_instances(staList);
