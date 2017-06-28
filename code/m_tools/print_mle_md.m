function print_mle_md(traceName,ctlg_m,ctlg_d,notes)

if (nargin<4)
	notes = '';
end

outDirName = 'fig/mle_md/';
details    = strcat(['m',num2str(ctlg_m),'_d',num2str(round(ctlg_d)),'km_',traceName,'_',notes,'.eps']);
fileName   = strcat([outDirName,details]);

set(gcf,'PaperPositionMode','auto') % --> ensures that sizes are the same as 
                                    % on screen
print('-depsc','-tiff','-r150',fileName)
