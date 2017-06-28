function print_traces(wformDir,magn,eD,eqName,nw,sta,chan)

h = figure(1);

name        = textscan(wformDir,'%s','delimiter','/');
dataSetName = name{1}{end};
outDirName  = strcat(['plots/',dataSetName,'/']);

if (~exist(outDirName,'dir'))
    ucmd = strcat(['mkdir ',outDirName]); unix(ucmd);
end

details     = strcat(['m',num2str(magn),'_d',num2str(eD),'_id',eqName, ...
    '_on_',nw,'.',sta,'.',chan,'.traces.eps']);
fileName = strcat([outDirName,details]);

set(gcf,'PaperPositionMode','auto') % --> ensures that sizes are the same as 
                                    % on screen
print('-depsc','-tiff','-r150',fileName)