function out = plot_ecdf_in_mBins(trList,valName,fig)
%Plot empirical cumulative distribution functions of quantity specified as
% <valName>.

global ftSize colours
m = trList.m;

if isfield(fig,'mRanges'); mRanges = fig.mRanges;
else                       mRanges = [6.5 8; 6 6.5; 5.5 6; 5 5.5; 4.5 5; 4 4.5; 3.5 4; 3 3.5; 2.5 3; 2 2.5];
end
mRanges  = mRanges(mRanges(:,1)>=min(trList.m),:);
nr       = size(mRanges,1);
%colours  = {'k','r','c','m','y','b',[0 .4 0],[255/256 165/256 0],[238/256 130/256 238/256],'k'};

hf=figure(fig.figNum); if fig.newFig; clf; end
hold on; grid on; box on; whitebg('w');

ph         = zeros(nr,1);
n          = zeros(nr,1);
F          = cell (nr,1);
X          = cell (nr,1);
Values     = cell(nr,1);
lgdStrings = cell(nr,1);
    
for ir = nr:-1:1; 

    fprintf(1,sprintf('%i/%i\n',ir,nr))
    
    mLo    = mRanges(ir,1);
    mUp    = mRanges(ir,2);
    idx    = find(m>=mLo & m<mUp);
    n(ir)  = numel(idx);
    
    if strcmp(valName,'hypDist')
        vals = trList.hypDist(idx);
        xlab = 'Hypocentral distance [km]';
    
    elseif strcmp(valName,'eqZ')
        vals = trList.eqZ(idx);
        xlab = 'Hypocenter depth [km]';
        
	elseif strcmp(valName,'noise')
        vals = cellfun(@(x) x(3),trList.accNoise(idx));
        xlab = 'Acceleration Noise [m/s/s]';
        
    else
        error(sprintf('Plotting ECDF for %s has not been implemented yet. Do it. Now.\n',valName))
        pause
    end
    
    % Compute & plot ecdf
    [f,x]  = ecdf(vals);
	ph(ir) = plot(x,f,'color',colours{ir},'lineWidth',2);

    % Save for output & legend
    F{ir}          = f;
    X{ir}          = x;
    Values{ir}     = vals;
    lgdStrings{ir} = sprintf('%3.1f <= M <%3.1f  (n=%i)',mRanges(ir,:),n(ir));;
end


%% Format Figure
%set(gca,'fontSize',ftSize,'YAxisLocation','right','ylim',[0 1])
set(gca,'fontSize',ftSize,'ylim',[0 1],'xscale',fig.yscale)
ylabel('\Phi (x \leq x'') ','fontSize',ftSize)
if isfield(fig,'valRange'); set(gca,'xlim',fig.valRange); end
xlabel(xlab,'fontSize',ftSize)


%% Legend
if fig.plotLegend; hl = legend(ph,lgdStrings);
                   set(hl,'location','NorthWest','fontSize',ftSize)
end

    
%% Output    
out.figHandle   = hf;
out.lineHandles = ph;
out.yEcdf       = F;
out.xEcdf       = X;
out.Values      = Values;