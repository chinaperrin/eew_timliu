% Plot attenuation curves for different magnitude bins, in order to find
% erroneous data.
%
% menandrin@gmail.com, 130624


%% 
addpath(genpath('../../../matlab/fct_downloads/'))
addpath(genpath('../../../../gmpe/'))
addpath(genpath('../../../../matlab/sac/'))
addpath(genpath('../../../../matlab/myFct/maps/'))
addpath(genpath('../../../../VS/data/plot/'))
clear all



%% ===============
%  0. Preparations

mmin = 3;
mmax = 9;
dm   = 0.5;
mm   = mmin:dm:mmax;
nm   = numel(mm)-1;

xband = 4;
xsnp  = 6;

alogmax = 5;
alogmin = -10;

level = 2;         % Level for percentiles

% Prepare distribution arrays
dint   = 20;
dedges = 0:dint:500;
nedges = numel(dedges)-1;
dd     = 0:500;



%%
% Station list
fprintf(1,'\n\n\n--------------------------------------------------------------------------------------')
fprintf(1,'\nImporting station lists\n')

% California
stFileName         = '../../../data/stations/eewvs-stations_allCali.hinv.fmtd2';
[stationList_cali] = import_stationlist_fB(stFileName);

stFileName     = '../../../data/stations/japan/sitepub_all_en.txt';
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
            
%dataSetNames2 = {'/scratch/memeier/VS/data/wform/CMS/scsn_M5to9_900101_011231/1band/'};

%dataSetNames3 = {'/scratch/memeier/VS/data/wform/socal/scsn_900101_011231/9bandsNew/'};
dataSetNames = {'/scratch/memeier/VS/data/wform/japan/k_kik/M6p_D25km/9bands/'};

[trL_sc,fc,snpLength,nsnp]    = import_fbOut(dataSetNames);
%[trL_jp,fc3,snpLength3,nsnp3] = import_fbOut(dataSetNames3);

fprintf(1,'Reading sensor orientations ... ')
orientation_socal = trL_sc.orntCode;
orientation_japan = trL_jp.orntCode;
fprintf(1,'done.\n')


%%
tL_idx_outUp = [];
tL_idx_outLo = [];
drel_outUp   = []; 
drel_outLo   = []; 
        
tL_idx_outUp_sc90 = [];
tL_idx_outLo_sc90 = [];
drel_outUp_sc90   = []; 
drel_outLo_sc90   = []; 


for im = 1:nm
    
    m_lo = mm(im);
    m_up = mm(im+1);
    
    idx_sc      = find( (trL_sc.m>=m_lo)      & (trL_sc.m<m_up) );
    %idx_sc90 = find( (trL_sc90.m>=m_lo) & (trL_sc90.m<m_up) );
    %idx_jp      = find( (trL_jp.m>=m_lo)      & (trL_jp.m<m_up) );
    
    
    % Prepare Plot
    figure(342); clf; hold on;
    grid on
    set(gca,'xscale','log')
    xlabel('Epicentral distance [km]','fontSize',15)
    ylabel('log_{10} max amplitudes [m/s]','fontSize',15)
    title(['Maximum velocities in ',num2str(snpLength),'sec since the p-pick, ' ...
        ,num2str(m_lo),' <= M < ',num2str(m_up)],'fontSize',15)
    set(gca,'xlim',[0.1 1000],'fontSize',15)
    set(gca,'ylim',[alogmin alogmax],'fontSize',15)
    set(gca,'ytick',alogmin:1:alogmax)
    whitebg(gcf,'k')
    set(gcf,'InvertHardcopy','off')
      
    
    % Socal   -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    if (~isempty(idx_sc))

        dVal = trL_sc.epiDist(idx_sc);
        aVal = trL_sc.pgv(idx_sc);
        ornt = orientation_socal(idx_sc);
        
%         dVal_N = dVal(strcmp(ornt,'N'));
%         aVal_N = aVal(strcmp(ornt,'N'));
%         h1=plot(dVal_N,log10(aVal_N),'xw');
% 
%         dVal_E = dVal(strcmp(ornt,'E'));
%         aVal_E = aVal(strcmp(ornt,'E'));
%         h2=plot(dVal_E,log10(aVal_E),'xc');
% 
%         dVal_Z = dVal(strcmp(ornt,'Z'));
%         aVal_Z = aVal(strcmp(ornt,'Z'));
%         h3=plot(dVal_Z,log10(aVal_Z),'xr');
        
        
        % Locate, plot and save outliers
        [~,up,~,~]   = CH2007_2sigma(m_up,dVal,'PGV','H','S','S',0);
        up           = up/1e2;
        diff         = aVal-up;
        idx_outUp    = find(diff>0);              
        tL_idx_outUp = [tL_idx_outUp; idx_sc(idx_outUp)];
        dtmp         = aVal(idx_outUp)./up(idx_outUp);    
        drel_outUp   = [drel_outUp; dtmp];

        
        [~,~,lo,~]   = CH2007_2sigma(m_lo,dVal,'PGV','H','S','R',0);
        lo           = lo/1e2;
        diff         = -aVal+lo;
        idx_outLo    = find(diff>0);              
        tL_idx_outLo = [tL_idx_outLo; idx_sc(idx_outLo)];
        dtmp         = aVal(idx_outLo)./lo(idx_outLo);    
        drel_outLo   = [drel_outLo; dtmp];
        
        h1 = plot(dVal,log10(aVal),'xw');
        
     end



%     % Socal 90ies   -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
%     if (~isempty(idx_sc90))
% 
%         dVal = trL_sc90.epiDist(idx_sc90);
%         aVal = cell2mat(trL_sc90.amax(idx_sc90));
%         ornt = orientation_socal90ies(idx_sc90);
% 
%         % Locate, plot and save outliers
%         [~,up,~,~]        = CH2007_2sigma(m_up,dVal,'PGV','H','S','S',0);
%         up                = up/1e2;
%         diff              = aVal-up;
%         idx_outUp         = find(diff>0);
%         tL_idx_outUp_sc90 = [tL_idx_outUp_sc90; idx_sc90(idx_outUp)];
%         dtmp              = aVal(idx_outUp)./up(idx_outUp);    
%         drel_outUp_sc90   = [drel_outUp_sc90; dtmp];
% 
%         
%         [~,~,lo,~]        = CH2007_2sigma(m_lo,dVal,'PGV','H','S','R',0);
%         lo                = lo/1e2;
%         diff              = aVal-lo;
%         idx_outLo         = find(diff<0);
%         tL_idx_outLo_sc90 = [tL_idx_outLo_sc90; idx_sc90(idx_outLo)];
%         dtmp              = aVal(idx_outLo)./lo(idx_outLo);    
%         drel_outLo_sc90   = [drel_outLo_sc90; dtmp];
%         
%           
%         h2 = plot(dVal,log10(aVal),'xy');
%         
% %         % Check if the indices with respect to <trL_sc90> are correct.
% %         hold on
% %         plot(dVal(idx_outLo),log10(aVal(idx_outLo)),'or')
% %         plot(trL_sc90.epiDist(idx_sc90(idx_outLo)),log10(cell2mat(trL_sc90.amax(idx_sc90(idx_outLo)))),'sb')
%         
%     end
%     
    
    
    
% 
%     % Japan   -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
%     if (~isempty(idx_jp))
% 
%         dVal = trL_jp.epiDist(idx_jp);
%         aVal = cell2mat(trL_jp.amax(idx_jp));
%         ornt = orientation_japan(idx_jp);
%         
%         dVal_N = dVal(strcmp(ornt,'N'));
%         aVal_N = aVal(strcmp(ornt,'N'));
%         h5=plot(dVal_N,log10(aVal_N),'ok');
% 
%         dVal_E = dVal(strcmp(ornt,'E'));
%         aVal_E = aVal(strcmp(ornt,'E'));
%         plot(dVal_E,log10(aVal_E),'oc')
% 
%         dVal_Z = dVal(strcmp(ornt,'Z'));
%         aVal_Z = aVal(strcmp(ornt,'Z'));
%         plot(dVal_Z,log10(aVal_Z),'or')
% 
%     end
%     
%     
%     if ( ~isempty(idx_sc) & ~isempty(idx_sc90) & ~isempty(idx_jp) )
%         legend([h1;h2;h3;h4;h5],'Socal North','Socal East','Socal Vertical','Socal90ies','Japan')    
%     elseif ( isempty(idx_sc) & isempty(idx_sc90) & isempty(idx_jp) )
%         flg_allEmpty=true;
%     end
%     
    %pause
    
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
              
%         legend([h1; h4; h4b; h4c],'data points, north component', ...
%             'CH2007, M = 2.75, rock',['CH2007, M=',num2str(m_up),', sigma_{up}, soil'], ...
%             ['CH2007, M=',num2str(m_lo),', sigma_{low}, rock'],'Location','NorthWest');
        print('-dpng','-r150',['gmpe_' num2str(im) '.png']);

        1+1;
        %pause
        
end



%% Add outliers to separate traceList and sort them

% SoCal - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (exist('idx_sc','var'))
    
    % Create sublist with outliers and sort such that largest errors come first
    OutUp                    = trL_sc.selectSubList(tL_idx_outUp);
    [drel_outUp_std,idx_std] = sort(drel_outUp,'descend');
    OutUp.sortList(idx_std);
    
    OutLo                    = trL_sc.selectSubList(tL_idx_outLo);
    [drel_outLo_std,idx_std] = sort(drel_outLo,'ascend');
    OutLo.sortList(idx_std);
    
    % Plot errors   
    figure(33);clf; hold on
    e1 = plot(drel_outLo_std,'ob');
    e2 = plot(drel_outUp_std,'xk');
    set(gca,'yscale','log')
    grid on
    l1 = legend([e1; e2],'Outliers << CH2007','Outliers >> CH2007','fonzSize',15);
    set(l1,'fontSize',15)
    ylabel('relative amplitude wrt CH2007 (aVal/CH07)','fontSize',15)
end



% SoCal90ies  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (exist('idx_sc90','var'))

    % Create sublist with outliers and sort such that largest errors come first
    OutUp_sc90                      = trL_sc90.selectSubList(tL_idx_outUp_sc90);
    [drel_outUp_sc90_std,sortIndex] = sort(drel_outUp_sc90,'descend');
    OutUp_sc90.sortList(sortIndex);
    
    OutLo_sc90                      = trL_sc90.selectSubList(tL_idx_outLo_sc90);
    [drel_outLo_sc90_std,sortIndex] = sort(drel_outLo_sc90,'ascend');
    OutLo_sc90.sortList(sortIndex);
    
    % Plot errors   - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    figure(34);clf; hold on
    e1 = plot(drel_outLo_sc90_std,'ob');
    e2 = plot(drel_outUp_sc90_std,'xk');
    set(gca,'yscale','log')
    grid on
    l1 = legend([e1; e2],'Outliers << CH2007','Outliers >> CH2007','fonzSize',15);
    set(l1,'fontSize',15)
    ylabel('relative amplitude wrt CH2007 (aVal/CH07)','fontSize',15)

end



% Japan       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




%% %%%%%%%%%%%%%%%%%%%%%%%%
%  Analyse oulier lists  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Preselection, e.g. to get rid of traces with only small errors
%OutLo_mod = OutLo.selectSubList(1:4000);


% Choose traceLists. Analyse all entries of <smallList> and use <largeList>
% to find e.g. co-recorded traces or traces from the same station.
largeList = trL_sc;
%largeList = trL_sc90;

%smallList = trL_sc90;
%smallList = OutLo_sc90;
%smallList = OutUp_sc90;
smallList = OutUp;
%smallList = OutLo_mod;


% A. Analyse station lists of outliers
% ====================================
n       = numel(smallList.sRate);
staList = cell(n,1);

for i = 1:n
    [~,~,~,sta_tmp,~] = get_stationCoords_cali(smallList.fullName{i},stationList_cali);
    staList{i}        = sta_tmp;
end
[countList] = count_station_instances(staList);

nsta = size(countList,1);
for ista = 1:nsta
   
    fprintf(1,['Looking at ',num2str(ista),'th of ',num2str(nsta),' stations \n'])
   
   stationName = countList{ista,1};
   [idx]        = find_station_traces(stationName,largeList);
   
   nout              = countList{ista,2};
   ntot              = numel(idx);
   countList{ista,3} = nout/ntot*100;       % Add fraction of outliers to countList
end



% B. Go through individual outlier traces
% =======================================

ntr = numel(smallList.sRate);
for itr = 1:ntr
    
    traceFullName = smallList.fullName{itr};
    slashIdx      = regexp(traceFullName,'/');
    recordName    = traceFullName(slashIdx(end)+1:end);
    fprintf(1,['\n\n\n     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n'])
    fprintf(1,['     Looking at ',recordName,'\n'])
    fprintf(1,['     . . . . . . . . . . . . . . . . .\n\n'])
    
    
    
    % Plot all amplitudes from the same station
    % -----------------------------------------
    fprintf(1,['Searching for traces from the same station ...\n'])
    
    % Make traceList with all traces of a station
    stationName              = recordName(end-10:end-8);
    [traceIndices]           = find_station_traces(stationName,largeList);
    trL_sameStation          = largeList.selectSubList(traceIndices);
    
    [ort,~] = find_orientation(trL_sameStation);
    lgc_z = cellfun(@(x) strcmp(x,'Z'), ort);
    lgc_n = cellfun(@(x) strcmp(x,'N'), ort);
    lgc_e = cellfun(@(x) strcmp(x,'E'), ort);
    
    figure(12); clf; hold on
    d = trL_sameStation.epiDist;
    a = log10(cell2mat(trL_sameStation.amax));
    m = trL_sameStation.m;
    scatter(d(lgc_z),a(lgc_z),50,m(lgc_z),'filled','marker','d')
    scatter(d(lgc_n),a(lgc_n),120,m(lgc_n),'filled','marker','p')
    scatter(d(lgc_e),a(lgc_e),50,m(lgc_e),'filled','marker','^')
    
    pd = plot(smallList.epiDist(itr),log10(smallList.amax{itr}), ...
        'pk','markerSize',14,'markerFaceColor','r');
    
    grid on
    colorbar
    title(['All traces of station ',stationName],'fontSize',15)
    set(gca,'xscale','log','xlim',[0.1 1000],'ylim',[alogmin alogmax],'ytick', ...
        alogmin:1:alogmax,'clim',[3 8])
    
    % GMPEs
    magnitude  = smallList.m(itr);
    [Ym,~,~,~] = CH2007(magnitude,dd,'PGV','H','S','R',0);
    pm         = plot(dd,log10(Ym/1e2),'-r','lineWidth',2);
    
    [Y3,~,~,~]  = CH2007(3,dd,'PGV','H','S','R',0);
    [Y4,~,~,~]  = CH2007(4,dd,'PGV','H','S','R',0);
    [Y5,~,~,~]  = CH2007(5,dd,'PGV','H','S','R',0);
    [Y6,~,~,~]  = CH2007(6,dd,'PGV','H','S','R',0);
    [Y7,~,~,~]  = CH2007(7,dd,'PGV','H','S','R',0);
    p3 = plot(dd,log10(Y3/1e2),'-k');
    p4 = plot(dd,log10(Y4/1e2),'--k');
    p5 = plot(dd,log10(Y5/1e2),':k');
    p6 = plot(dd,log10(Y6/1e2),'-b');
    p7 = plot(dd,log10(Y7/1e2),'--b');
    l1 = legend([pd;pm;p7;p6;p5;p4;p3],'Current data point',['CH2007 (M',num2str(magnitude),')'],'M7','M6','M5','M4','M3');
    set(l1,'location','SouthWest')
    
    
    
    % Plot all amplitudes of same event
    % ---------------------------------
    fprintf(1,['Searching for traces from the same event ...\n'])

    ptIdx         = regexp(recordName,'\.');
    eventName     = recordName(1:ptIdx(1)-1);

    [traceIndices] = find_event_traces(eventName,largeList);
    trL_sameEvent  = largeList.selectSubList(traceIndices);
    
    figure(13); clf; hold on
    d = trL_sameEvent.epiDist;
    a = log10(cell2mat(trL_sameEvent.amax));
    m = trL_sameEvent.m;
    scatter(d,a,30,m,'filled')
    pd = plot(smallList.epiDist(itr),log10(smallList.amax{itr}), ...
        'pk','markerSize',14,'markerFaceColor','r');
    
    grid on
    colorbar
    title(['All traces of event ',eventName],'fontSize',15)
    set(gca,'xscale','log','xlim',[0.1 1000],'ylim',[alogmin alogmax],'ytick', ...
        alogmin:1:alogmax,'clim',[3 8])
    
    % GMPEs
    magnitude  = smallList.m(itr);
    [Ym,~,~,~] = CH2007(magnitude,dd,'PGV','H','S','R',0);
    pm         = plot(dd,log10(Ym/1e2),'-r','lineWidth',2);
    
    [Y3,~,~,~]  = CH2007(3,dd,'PGV','H','S','R',0);
    [Y4,~,~,~]  = CH2007(4,dd,'PGV','H','S','R',0);
    [Y5,~,~,~]  = CH2007(5,dd,'PGV','H','S','R',0);
    [Y6,~,~,~]  = CH2007(6,dd,'PGV','H','S','R',0);
    [Y7,~,~,~]  = CH2007(7,dd,'PGV','H','S','R',0);
    p3 = plot(dd,log10(Y3/1e2),'-k');
    p4 = plot(dd,log10(Y4/1e2),'--k');
    p5 = plot(dd,log10(Y5/1e2),':k');
    p6 = plot(dd,log10(Y6/1e2),'-b');
    p7 = plot(dd,log10(Y7/1e2),'--b');
    l2 = legend([pd;pm;p7;p6;p5;p4;p3],'Current data point',['CH2007 (M',num2str(magnitude),')'],'M7','M6','M5','M4','M3');
    set(l2,'location','SouthWest')

    
    % Plot waveforms
    % ------------------------
    plotTraceSummary(traceFullName,largeList,60,1,10,1)

%    if (strcmp(recordName(end-4),'Z'))
        1+1;
%    end
end





%%



%n       = 3000;
n       = numel(OutLo.sRate);
staList = cell(n,1);

for i = 1:n
    %[~,~,~,sta_tmp,~] = get_stationCoords_cali(OutUp.fullName{i},stationList_cali);
    [~,~,~,sta_tmp,~] = get_stationCoords_cali(OutLo.fullName{i},stationList_cali);
    staList{i}        = sta_tmp;
end
[countList] = count_station_instances(staList);



nsta = 10;
for i = 1:nsta
    
    % Make traceList with all traces of a station
    [traceIndices]           = find_station_traces(countList{i},trL_sc);
    flg_select               = false(numel(trL_sc.t0),1);
    flg_select(traceIndices) = true;
    trL_tmp                  = trL_sc.selectSubList(flg_select);
    
    figure(12); clf; hold on
    d = trL_tmp.epiDist;
    a = log10(cell2mat(trL_tmp.amax));
    m = trL_tmp.m;
    scatter(d,a,10,m,'filled')
    set(gca,'xscale','log')
    set(gca,'xlim',[0.1 1000])
    %set(gca,'ylim',[alogmin alogmax])
    set(gca,'ytick',[alogmin:1:alogmax])
    grid on
    colorbar
    title(['Station ',countList{i,1},' (',num2str(countList{i,2}),' outliers)'],'fontSize',15)
    
    % GMPEs
    [Y3,~,~,~]  = CH2007(3,dd,'PGV','H','S','R',0);
    [Y4,~,~,~]  = CH2007(4,dd,'PGV','H','S','R',0);
    [Y5,~,~,~]  = CH2007(5,dd,'PGV','H','S','R',0);
    [Y6,~,~,~]  = CH2007(6,dd,'PGV','H','S','R',0);
    [Y7,~,~,~]  = CH2007(7,dd,'PGV','H','S','R',0);
    p3=plot(dd,log10(Y3/1e2),'-k');
    p4=plot(dd,log10(Y4/1e2),'--k');
    p5=plot(dd,log10(Y5/1e2),':k');
    p6=plot(dd,log10(Y6/1e2),'-b');
    p7=plot(dd,log10(Y7/1e2),'--b');
    legend([p3;p4;p5;p6;p7],'M3','M4','M5','M6','M7')
    pause
end


idx=find(cell2mat(trL_tmp.amax)<5e-7);
plot(trL_tmp.epiDist(idx),log10(cell2mat(trL_tmp.amax(idx))),'xr');


% Traces with amax = 0  - - - - - - - - - - - - - - - - - - - - - - - - - -  
idx_zilt   = find(cell2mat(OutLo.amax)==0);
OutLo_zilt = OutLo.selectSubList(idx_zilt);

n       = numel(OutLo_zilt.t0);
staList = cell(n,1);
for i = 1:n
    [~,~,~,sta_tmp,~] = get_stationCoords_cali(OutLo_zilt.fullName{i},stationList_cali);
    staList{i}        = sta_tmp;
end
[countList] = count_station_instances(staList);
