function [hf] = cfRegressionModels(regObjects,statistic)

global fc 

ftSize = 12;

nm      = numel(regObjects);    % No. of models to compare
nbands  = size(fc,1);
xb      = 1:nbands;
Symbols = {'o','x','^','v','d','s','p','*'};
Colors  = {'y','w','c','r','g','b'};

if strcmp(statistic,'bic');    figNr = 121; end
if strcmp(statistic,'sigma2'); figNr = 122; end

hf = figure(figNr); clf;
ha = tight_subplot(1,2,[0.01 0.05],[0.1 0.05],[0.12 0.1]);
whitebg(gcf,'k')
set(gcf,'InvertHardcopy','off')

legendString = cell(nm,1);
H            = [];
for im = 1:nm
    
    REG = regObjects{im};
    if strcmp(statistic,'bic')
       zStat = cell2mat(REG.zBic);
       hStat = cell2mat(REG.hBic);
    elseif strcmp(statistic,'sigma2')
       zStat = cell2mat(REG.zSigma2);
       hStat = cell2mat(REG.hSigma2);
    end
    
    symbol = Symbols{im};
    color  = Colors{im};
    
    % Vertical
    axes(ha(1)); grid on; box on; hold on
    hz = plot(xb,zStat,'marker',symbol,'lineStyle','none','color',color,'lineWidth',2);
    H  = [H,hz];
    ylabel(statistic,'fontSize',ftSize)
    if (im==1); title([statistic,' on Z-comps'],'fontSize',ftSize);
                xlabel('Frequency band','fontSize',ftSize); 
                set(gca,'xTick',xb,'xTickLabel',xb,'fontSize',ftSize);
                set(gca,'fontSize',ftSize); end 
    
    % Horizontal
    axes(ha(2)); grid on; box on; hold on
    hz = plot(xb,hStat,'marker',symbol,'lineStyle','none','color',color,'lineWidth',2);
    if (im==1); title([statistic,' on H-comps'],'fontSize',ftSize);
                xlabel('Frequency band','fontSize',ftSize); 
                set(gca,'xTick',xb,'xTickLabel',xb,'fontSize',ftSize);
                set(gca,'fontSize',ftSize); end 
    
    legendString{im} = REG.modelName;
end

l1 = legend(H,legendString,'Location','North');
set(l1,'fontSize',ftSize)