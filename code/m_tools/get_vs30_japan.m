function [vs30] = get_vs30_japan(stationList,opt)
% READ VS30 VALUE FOR LIST OF STATIONS FROM STATION-FILES
% For k- & kik-Net there is a vs and vp profile available for each site.
% They are stored in eq/data/japan/site/kik/ & .../knt/

vs30_default = 500;
zMaster      = (1:30)';
ftSize       = 15;

% Read all sac-files in benchmarking directory
kikDir       = '~/Documents/eq/data/japan/site/kik/';
fileList    = dir(kikDir);                % get a list of all files
fileNameList = {fileList.name}';
datIdx      = regexp(fileNameList,'\.dat$');
isDat       = cellfun(@(x) ~isempty(x), datIdx);
kikFileList = fileNameList(isDat);

kntDir       = '~/Documents/eq/data/japan/site/knt/';
fileList    = dir(kntDir);                % get a list of all files
fileNameList = {fileList.name}';
datIdx      = regexp(fileNameList,'\.dat$');
isDat       = cellfun(@(x) ~isempty(x), datIdx);
kntFileList = fileNameList(isDat);

% Check if some station names exist in both files? No.
% ct = 0;
% for ist = 1:numel(kntFileList)
%     thisName = kntFileList{ist};
%     idxHits  = find(cellfun(@(x) ~isempty(x),  regexp(thisName,kikFileList)));
%     if ~isempty(idxHits); ct=ct+1; end    
% end

nst             = numel(stationList.lat);
vs30.meanSlow   = zeros(nst,1);
vs30.medianSlow = zeros(nst,1);
vs30.meanVs     = zeros(nst,1);
vs30.maxDepth   = zeros(nst,1);
ct0 = 0;
ct2 = 0;
for ist = 1:nst

    print_iteration_numbers(ist,nst,'hundreds')
    
    % Find station 
    stName      = stationList.name{ist};
    idxKikHits  = find(cellfun(@(x) ~isempty(x), regexp(kikFileList,stName)));
    idxKntHits  = find(cellfun(@(x) ~isempty(x), regexp(kntFileList,stName)));
    
    nhits       = numel(idxKikHits)+numel(idxKntHits);
    if nhits==1
        
        % Find file
        if     ~isempty(idxKikHits) & isempty(idxKntHits); fileFullName = sprintf('%s%s',kikDir,kikFileList{idxKikHits});
        elseif  isempty(idxKikHits) &~isempty(idxKntHits); fileFullName = sprintf('%s%s',kntDir,kntFileList{idxKntHits});
        end
        
        % Read profile files
        vsProfile = read_bosai_depth_profile(fileFullName,opt);
    
        % Check if need to use default values
        if vsProfile.fileType==0 |unique(vsProfile.vs)==0
            vsProfile.vs = vs30_default*ones(size(zMaster));
            vsProfile.z  = zMaster;
        end

        if max(vsProfile.z)<30;
            if vsProfile.vs(end)<vs30_default
                vsProfile.vs = [vsProfile.vs; vs30_default];
            else
                vsProfile.vs = [vsProfile.vs; vsProfile.vs(end)];
            end
            vsProfile.vp = [vsProfile.vp; vsProfile.vs(end)*1.7];
            vsProfile.z  = [vsProfile.z ; 30];
        end
        
        % Interpolate slowness
        slowness     = 1./vsProfile.vs;
        slowness(slowness==inf) = nan;
        %slowness_itp = interp1(vsProfile.z,slowness,zMaster,'spline');
        slowness_itp = interp1(vsProfile.z,slowness,zMaster,'linear');

        vs30.meanSlow  (ist) = 1/mean(slowness_itp);
        vs30.medianSlow(ist) = 1/median(slowness_itp);
        vs30.meanVs    (ist) = mean(vsProfile.vs);
        
        vs30.maxDepth(ist)   = max(vsProfile.z);
        
        if max(vsProfile.vs)>1500; xmax=max(vsProfile.vs);
        else                    xmax=1500;
        end
        slashIdx = regexp(fileFullName,'/');
        tString  = fileFullName(slashIdx(end-3)+1:end);
        
        if opt.plotVsProfile
        %if opt.plotVsProfile &vsProfile.fileType==2
            clf; hold on; grid on; box on;
            p1=plot(1./slowness,vsProfile.z,'-or','lineWidth',1,'markerSize',9); 
            p2=plot(1./slowness_itp,zMaster,'-xk'); 
            set(gca,'ydir','reverse','xlim',[0 xmax],'ylim',[0 35]) 
            l1=line([vs30.medianSlow(ist) vs30.medianSlow(ist)],get(gca,'ylim'),'lineWidth',1)
            l2=line([vs30.meanSlow(ist)   vs30.meanSlow(ist)  ],get(gca,'ylim'),'lineWidth',2,'lineStyle',':' ,'color','r')
            l3=line([vs30.meanVs(ist)     vs30.meanVs(ist)    ],get(gca,'ylim'),'lineStyle','-.','color','m')
            lg1 = legend([p1;p2;l1;l2;l3],'Tabulated profile + last value','Interpolated profile','Vs30 (median slowness)','Vs30 (mean slowness)','Vs30 (mean Vs)');
            title(tString)
            xlabel('Vs [m/s]','fontSize',ftSize,'fontName','Avenir')
            ylabel('Depth [m]','fontSize',ftSize,'fontName','Avenir')
            set(gca,'fontSize',ftSize,'fontName','Avenir') 
            if opt.printVsProfile; print('-dpng',strrep(fileFullName,'.dat','.png')); end
        end
        

    else
        if isempty(idxKikHits) & isempty(idxKntHits)
            fprintf(1,sprintf('No depth profile file found for %s\n',stName))
            ct0 = ct0+1;
            vs30.meanSlow  (ist) = vs30_default;
            vs30.medianSlow(ist) = vs30_default;
            vs30.meanVs    (ist) = vs30_default;
            
        elseif ~isempty(idxKikHits) &~isempty(idxKntHits)
            fprintf(1,sprintf('2 depth profile file found for %s\n',stName))
            ct2 = ct2+1;
        end
    end
end
fprintf(1,sprintf('\nNo depth profiles found in %i cases\n',ct0))
fprintf(1,sprintf('Depth profiles found in both k- & kik-List in %i cases\n',ct2))


if opt.plotVsMap
    stLatVect = stationList.lat;
    stLonVect = stationList.lon;
    vs30Vect  = vs30.meanSlow;
    
    figure(2119); clf; hold on; grid on; box on;
    scatter(stLonVect,stLatVect,15,vs30Vect,'filled')
    set(gca,'xlim',[123 147],'ylim',[23 46])
    colormap(hot)
    colorbar
    caxis([180 760])
    % print('-dpng','~/Documents/eq/data/japan/site/map/bosai_vs30map.png')
    
%     addpath(genpath('../../../matlab/colorscales/'))
%     clrIncs  = 170:10:800;
%     clrSteps = [0;180;360;760;2000];
%     clrSteps = [0;180;240;300;360;490;620;760;2000];
%     clrmap   = make_cScale_w_irregular_steps(clrIncs,clrSteps);
%     colorbar; colormap(clrmap)
%     colormap(hot)
end