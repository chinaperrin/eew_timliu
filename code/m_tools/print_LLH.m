function print_LLH(magn,eD,eqName,nw,sta,chan)

h = figure(333);

outDirName  = strcat(['plots/mle/']);

if (~exist(outDirName,'dir'))
    ucmd = strcat(['mkdir ',outDirName]); unix(ucmd);
end

details     = strcat(['m',num2str(magn),'_d',num2str(eD),'_id',eqName, ...
    '_on_',nw,'.',sta,'.',chan,'.LLH.eps']);
fileName = strcat([outDirName,details]);

set(gcf,'PaperPositionMode','auto') % --> ensures that sizes are the same as 
                                    % on screen
print('-depsc','-tiff','-r150',fileName)