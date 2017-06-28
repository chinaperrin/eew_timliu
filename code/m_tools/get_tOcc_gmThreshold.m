function [tOcc] = get_tOcc_gmThreshold(trList,gmMeasure,itFieldList)

nt     = numel(itFieldList);
tOcc.t = cell(nt,1);
tOcc.m = cell(nt,1);
tOcc.r = cell(nt,1);

for it = 1:nt

    itField    = itFieldList(it);
    
    if     strcmp(gmMeasure,'acc'); idxThd = cellfun(@(x) x.acc_zen (itField), trList.scalFeature);
    elseif strcmp(gmMeasure,'vel'); idxThd = cellfun(@(x) x.vel_zen (itField), trList.scalFeature);
    elseif strcmp(gmMeasure,'mmi'); idxThd = cellfun(@(x) x.immi_zen(itField), trList.scalFeature);
    elseif strcmp(gmMeasure,'pgx');
        if     it==1;               idxThd = trList.ppga.idx_zen_snp;
        elseif it==2;               idxThd = trList.ppgv.idx_zen_snp;
        elseif it==3;               idxThd = trList.ppgd.idx_zen_snp;
        end
    else   error('Invalid choice of gmMeasure.')
    end
    
    if ~isempty(idxThd)
        idx        = find(idxThd~=0);
        idxThd     = idxThd(idx);
        tmpList    = trList.selectSubList(idx);
        sr         = tmpList.sRate;
        tOcc.t{it} = idxThd./sr;
        tOcc.r{it} = tmpList.hypDist;
        tOcc.m{it} = tmpList.m;
    end
end






% % 0.01g
% has0p01g = cell2mat(cellfun(@(x) ~isempty(x.i0p01g), trList.scalFeature,'uniformOutput',0));
% tmpList  = trList.selectSubList(find(has0p01g));
% i0p01g   = cell2mat(cellfun(@(x) x.i0p01g, tmpList.scalFeature,'uniformOutput',0));
% sr       = tmpList.sRate;
% t0p01g   = i0p01g./sr;
% r        = tmpList.hypDist;
% p1 = plot(r,t0p01g,'ok','markerSize',mkSize,'markerFaceColor','k'); set(gca,'ylim',[0 20]);
% out.acc0p01g.t    = t0p01g;
% out.acc0p01g.r    = r;
% out.acc0p01g.m    = tmpList.m;
% out.acc0p01g.tppx = tmpList.tppx;
% 
% % 0.1g
% has0p1g = cell2mat(cellfun(@(x) ~isempty(x.i0p1g), trList.scalFeature,'uniformOutput',0));
% tmpList  = trList.selectSubList(find(has0p1g));
% i0p1g    = cell2mat(cellfun(@(x) x.i0p1g, tmpList.scalFeature,'uniformOutput',0));
% sr       = tmpList.sRate;
% t0p1g   = i0p1g./sr;
% r        = tmpList.hypDist;
% p2 = plot(r,t0p1g,'ok','markerSize',mkSize+1,'markerFaceColor','b'); set(gca,'ylim',[0 20]);
% out.acc0p1g.t = t0p1g;
% out.acc0p1g.r = r;
% out.acc0p1g.m = tmpList.m;
% out.acc0p1g.tppx = tmpList.tppx;
% 
% 
% % 0.5g
% has0p5g  = cell2mat(cellfun(@(x) ~isempty(x.i0p5g) , trList.scalFeature,'uniformOutput',0));
% tmpList  = trList.selectSubList(find(has0p5g));
% i0p5g    = cell2mat(cellfun(@(x) x.i0p5g, tmpList.scalFeature,'uniformOutput',0));
% sr       = tmpList.sRate;
% t0p5g   = i0p5g./sr;
% r        = tmpList.hypDist;
% p3 = plot(r,t0p5g,'ok','markerSize',mkSize+3,'markerFaceColor','c'); set(gca,'ylim',[0 20]);
% out.acc0p5g.t = t0p5g;
% out.acc0p5g.r = r;
% out.acc0p5g.m = tmpList.m;
% out.acc0p5g.tppx = tmpList.tppx;
% 
% % 1.0g
% has1p0g  = cell2mat(cellfun(@(x) ~isempty(x.i1p0g) , trList.scalFeature,'uniformOutput',0));
% tmpList  = trList.selectSubList(find(has1p0g));
% i1p0g    = cell2mat(cellfun(@(x) x.i1p0g, tmpList.scalFeature,'uniformOutput',0));
% sr       = tmpList.sRate;
% t1p0g   = i1p0g./sr;
% r        = tmpList.hypDist;
% p4 = plot(r,t1p0g,'ok','markerSize',mkSize+5,'markerFaceColor','r'); set(gca,'ylim',[0 20]);
% out.acc1p0g.t = t1p0g;
% out.acc1p0g.r = r;
% out.acc1p0g.m = tmpList.m;
% out.acc1p0g.tppx = tmpList.tppx;

% l1 = legend([p1; p2; p3; p4],'0.01g','0.1g','0.5g','1.0g');
% set(l1,'fontSize',ftSize,'fontName','Avenir')
% set(gca,'xlim',[0 50],'fontSize',ftSize,'fontName','Avenir')
% xlabel('Hypocentral Distance [km]','fontSize',ftSize,'fontName','Avenir')
% ylabel('Time since P-onset [s]','fontSize',ftSize,'fontName','Avenir')
% title('Occurence time of threshold acceleration [s]','fontSize',ftSize,'fontName','Avenir')

%print('-dpng','~/programs/filterBank/fig/i37/tpgx/new/tOcc_acceleration')