function analyse_sublist_wforms(trList)
%% Options und so ...
addpath(genpath('../m_tools/'))
addpath(genpath('../'))

o.scp_wforms = 1;
snpLength    = 0.5;


% Picker parameters
px.pxThreshold = 5;
px.fLow_px     = 3;
px.fUp_px      = 10;
px.fOrder      = 4;
px.ntap        = 100;            % Don't pick before ntap
px.sWindow     = 1;              % Window lenght for measuring signal amps after p-pick [sec]
px.nWindow     = 1;              % Window lenght for measuring noise before p-pick [sec]
px.wmax        = 1e1;
px.wmin        = .1;
px.sWtFct      = 'exp';          % 'lin' or 'exp' or '1-exp'
px.nWtFct      = 'exp';          % 'lin' or 'exp' or '1-exp'
px.plotPx      = true;
px.animatePx   = true;


%% Go through all traces
ntr = size(trList.fullName,1);


for itr = 1:ntr

    fullName                = trList.fullName{itr};
    [pathName,recordName,~] = splitFullName(fullName);
    threeCompPattern        = strcat([pathName,recordName,'*']);

    fprintf(1,sprintf('\n%i/%i: \trecord %s\n',itr,ntr,recordName))
    fprintf(1,['\t',num2str(trList.m(itr)),'M @ ',num2str(trList.hypDist(itr)),'km\n'])
    fprintf(1,sprintf('\tSNR: %6.2f\n ',trList.snr(itr)))
    
    if ( (~exist(trList.fullName{itr},'file')) && (o.scp_wforms) )
        scp_wform(threeCompPattern,pathName)
    end
    
    
    
    plot_traceSummary(fullName,trList,15,1,10,0) 
    1+1;
    
    % Manually call repick function
    % [ppxIdx2,snr2] = repick_SBPx(fullName,trList,tppx0,px)
    % trList.ppxIdx  = ppxIdx2;
    % trList.snr     = snr2;
    %
    % add2pxList(fullName,ppxIdx2,snr2,'','pxList_i34.mat')
    
end












%% APPENDIX
    % Plot various figures
    %[h1,h2,h3,h4] = reproduce_fbOut_1trace_minimal(i,tmpList,snpLength);
    
% 
%     tppx1 = tmpList.tppx(i)-dt1(i);     % Pick from simple picker
%     tppx2 = tmpList.tppx(i)-dt2(i);     % Pick from searchOpt picker
%     
%     
%     figure(h1)
%     hold on
%     ymax = abs(get(gca,'ylim')); 
%     ymax = ymax(1);
%     l1 = line([tppx1 tppx1],[-ymax ymax],'color','k','lineStyle','--');
%     l2 = line([tppx2 tppx2],[-ymax ymax],'color','r','lineStyle','-.');
%     title(['simplePx-SNR: ',num2str(snr1(i),'%4.0f'),'   ---   searchOptPx-SNR: ',num2str(snr2(i),'%4.0f')],'fontSize',ftSize)
%     
%     if (status==0)
%         % File found
%         fprintf(1,[num2str(i),'. Waveform-file found\n'])
%         fprintf(1,['   ',num2str(tmpList.m(i)),'M @ ',num2str(tmpList.hypDist(i)),'km\n'])
%         [h1,h2] = reproduce_fbOut_1trace_minimal(tmpList.fullName{i},tmpList,0.5,20,tmpList.orntCode{i});
%     else
%         fprintf(1,[num2str(i),'. Waveform-file has already been blocked (but amplitudes are still in saved traceLists)\n'])
%     end


% Low amp traces in Japan - - - - - - - - - - - - - - - - - - - - - - - - -
% atmp    = cell2mat(cellfun(@(x) x(5,snippet), TraceList.amax,'UniformOutput',0));
% idx     = find( (TraceList.m>5.5) & (log10(atmp)<-4) );
% tmpList = TraceList.selectSubList(idx);
% plot_data(tmpList.hypDist,tmpList.m,tmpList.amax,snippet,order);
% 
% figure; hold on;
% latmin=30; latmax=45; lonmin=128; lonmax=145;
% [p1] = plot_japan_map(latmin,latmax,lonmin,lonmax);
% grid on
% 
% p2 = plot(tmpList.stationLon,tmpList.stationLat,'ob','lineWidth',2);
% p3 = plot(tmpList.eqLon,tmpList.eqLat,'xr','lineWidth',2);
% l1=legend([p1,p2,p3],'coastline','station','hypocenter');
% set(l1,'fontSize',16)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


% Compare picks
% dt1  = cell2mat(tmpList.var);
% snr1 = cell2mat(tmpList.var3);
% 
% dt2  = cell2mat(tmpList.var2);
% snr2 = cell2mat(tmpList.var4);
