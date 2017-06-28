function [deltasVect,maxNoise] = get_median_pick_delay(acc,meta,n,k,nrnd,px,opts)

global fOrder ftSize

% Unpack parameters
ntap     = px.Param.ntap;
fLow_px  = px.Param.fLow_px;
fUp_px   = px.Param.fUp_px;
tw1      = meta.sr/px.Param.tw1_denom;
%tw2      = meta.sr/px.Param.tw2_denom;
 

% Z: Read & process VERTICAL trace ............................
t             = meta.t;
ns            = numel(acc);
ppxIdx        = meta.ppxIdx;
sr            = meta.sr;
dt            = 1/sr;


% Determine interval <idx0Vect> from which onset indices can be sampled:
% OUTDATED GRAPH
%
%                                                tPrime                 tppx
% | -------- | -----------------------------| ------------ | ------------ | 
%    ntap                 idxOnset                ns2           ns_gap    
%
% 1        ntap                           ppxIdx         ppxIdx         ppxIdx
%                                        - ns_gap       - ns_gap
%                                        - ns2

ns_gap   = ceil(sr/4); % No. of samples between noise window and pick that are not used
aNoise   = acc(1:ppxIdx-ns_gap);
maxNoise = prctile(aNoise(ntap:end),95);
tt       = (0:dt:.5)';
nt       = numel(tt);
idxOnset = ntap+tw1:1:ppxIdx-ns_gap-nt;    % Possible start points for synthetic signal:
                                       % After ntap, but before noise gap -
if ~isempty(idxOnset)
    
    idx0Vect   = randi([idxOnset(1),idxOnset(end)],nrnd,1);   % not after ns-gap before true onset.
    deltasVect = zeros(nrnd,1);
    
    for irnd =  1:nrnd
        
        idx0 = idx0Vect(irnd);
        yy   = [zeros(idx0,1); k*tt.^n;];  % From 1:idx0, amp of synthetic signal = 0
        
        % Make yy have length ns
        if length(yy)>ns; yy       = yy(1:ns); end
        if length(yy)<ns; nmissing = ns - length(yy);
                          yy       = [yy; yy(end)*ones(nmissing,1)];
        end
        
        % Make aNoise have length ns
        if length(aNoise)>ns; aNoise   = aNoise(1:ns); end
        if length(aNoise)<ns; nmissing = ns - length(aNoise);
                              aNoise   = [aNoise; aNoise(end)*ones(nmissing,1)];
        end

        y = aNoise+yy;
 
        % Filter and pick synthetic signal
        yb  = bworth(y,sr,[fLow_px,fUp_px],'band',fOrder,'causal');       % Bandpass filter for picking
        %         yb2 = bworth(y,sr,[1,fUp_px],'band',fOrder,'causal');       % Bandpass filter for picking
        %         yb3 = bworth(y,sr,[3,fUp_px],'band',fOrder,'causal');       % Bandpass filter for picking
        %         clf; hold on; plot(y); plot(yb,'r');  plot(yb2,'m');  plot(yb3,'k');
        
        % Modify filtered synthetic
        [maxVal,maxIdx] = max(yb);
        yb(maxIdx:end)  = maxVal;
        %plot(yb,'g')
        
        ppxIdx0 = idx0-100;
        if ppxIdx0+tw1<ntap; ppxIdx0=ntap+tw1; end;
        [ppxIdxS,~,~,~] = SBPx(yb,1/sr,ppxIdx0,px.Param,px.Weight,px.Opt);
        
        if ~isempty(ppxIdxS) 
            
            deltasVect(irnd) = ppxIdxS-idx0;
             
            if opts.plotSynth
                figure(4); clf; hold on; grid on; box on;
                plot(y,'r','lineWidth',2)
                plot(idx0,0,'ok','lineWidth',2,'markerSize',12,'markerFaceColor','w')
                plot(idx0,0,'or','lineWidth',2,'markerSize',10,'markerFaceColor','y')
                xlm = [1 idx0+200];
                line([xlm(1) xlm(2)], [maxNoise, maxNoise],'color','b','lineWidth',2,'lineStyle','-.');
                set(gca,'xlim',xlm,'ylim',[-2*maxNoise 2*maxNoise],'fontSize',ftSize)
                plot(yb,'-.b')
                
                plot(ppxIdxS,0,'dk','lineWidth',4,'markerSize',12)
                plot(ppxIdxS,0,'dw','lineWidth',1,'markerSize',12)
                %plot(t(ppxIdxS),0,'dy','lineWidth',2,'markerSize',12)
                %set(gca,'xlim',[ppxIdxS-50 ppxIdxS+100],'ylim',[-3e-3 3e-3],'fontSize',ftSize)
                set(gca,'xlim',[ppxIdxS-100 ppxIdxS+150],'ylim',[-3e-2 3e-2],'fontSize',ftSize)
                1+1;
                
                % Plot figure to explain pick-delay correction
                neverTrue=false;
                if neverTrue
                    figure(449); clf; hold on; grid on; box on;
                    mkSize = 18; 
                    yplot = k*tt.^n;
                    tplot = t(ppxIdx)-5.05+tt;
                    c1 = plot(t ,y ,'r','lineWidth',2);
                    c2 = plot(tplot,yplot       ,'--k','lineWidth',2);
                         plot([0 tplot(1)],[0 0],'--k' ,'lineWidth',2)
                    p2 = plot(tplot(1),0  ,'pk','markerFaceColor','y'     ,'markerSize',mkSize)
                    p1 = plot(8.9,1.5e-4,'pk','markerFaceColor',[0 .4 0],'markerSize',mkSize)
                    set(gca,'xlim',[8.3 9.1],'ylim',[-1e-4 5e-4],'fontSize',ftSize)
                    xlabel('Time since origin time [sec]','fontSize',ftSize)
                    ylabel('Acceleration [m/s/s]','fontSize',ftSize)
                    
                    %line([10.96 10.96],[-1.5e-4 1.5e-4])
                    l1 = legend([c1; c2; p1; p2],'Seismogram','y = t^4','original pick','corrected pick');
                    set(l1,'location','northWest','fontSize',ftSize)
                    
                    outFullName = '~/programs/filterBank/fig/i35/pickCorr/new/pickCorr_example.eps';
                    print('-depsc2',outFullName)
                end
            end 
        else
            deltasVect(irnd) = 99999; 
        end
        1+1;
    end
else
    deltasVect = 99999;
end



%% APPENDIX



% tPrime   = maxNoise^(1/n);    % Time for y=kt**n to reach maximum noise amp 
% if tPrime<0.02; tPrime=0.02; end
% tt       = (0:dt:10*tPrime)';
% yy       = k*tt.^n;
% ny       = numel(tt);         % No. of symples of entire synthetic function
% 
% idxOnset = ntap+tw1:1:ppxIdx-ns_gap-ny;    % Possible start points for synthetic signal:
%                                        % After ntap, but before noise gap -
% if ~isempty(idxOnset)
%     
%     idx0Vect   = randi([idxOnset(1),idxOnset(end)],nrnd,1);   % not after ns-gap before true onset.
%     deltasVect = zeros(nrnd,1);
%     
%     for irnd =  1:nrnd
%         
%         idx0 = idx0Vect(irnd);
%         y    = zeros(ns,1);
%         
%         % Interval 1: Only noise
%         sIdx1 = 1;
%         eIdx1 = idx0-1;
%         int1  = (sIdx1:eIdx1)';
%    
%         % Interval 2: Signal + noise
%         sIdx2 = idx0;               
%         eIdx2 = idx0+ny-1;
%         int2  = (sIdx2:eIdx2)';
% 
%         % Interval 3: Zeros
%         sIdx3 = eIdx2+1;
%         eIdx3 = ns;
%         int3  = (sIdx3:eIdx3)';
%         
%         %interval = [int1; int2; int3;];
%         y(int1) = aNoise(int1);
%         y(int2) = aNoise(int2)+yy;
%         y(int3) = 0;
% 
% 
%         
%         % Create snthetic signal
%         sIdx1 = 1;                  % Interval 1: Only noise
%         eIdx1 = idx0-1;
%         int1  = (sIdx1:eIdx1)';
%         
%         sIdx2 = idx0;               % Interval 2: Signal + noise
%         eIdx2 = idx0+ns2-1;
%         int2  = (sIdx2:eIdx2)';
%         
%         sIdx3 = eIdx2+1;            % Interval 3: Only signal
%         eIdx3 = sIdx3+ns3;
%         int3  = (sIdx3:eIdx3)';
%         
%         int4  = (eIdx3+1:eIdx3+1000)'; % Interval 4: Zeros
%         
%         y(int1)  = aNoise(int1);
%         y(int2)  = aNoise(int2)+yy(1:ns2);
%         y(int3)  = yy(ns2+1:ns2+ns3);
%         y(int4)  = 0;
%         
%         interval = [int1; int2; int3; int4];
%         if interval(end)>ns; interval = int1(1):1:ns; end
%         tn       = t(interval);