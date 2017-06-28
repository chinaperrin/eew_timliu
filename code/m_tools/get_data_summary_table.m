function [datTable] = get_data_summary_table(TraceList)

idx_scsn = find( strcmp(TraceList.dataSetName,'scsn') |strcmp(TraceList.dataSetName,'scsnPx') );
idx_jp   = find( strcmp(TraceList.dataSetName,'kNet') |strcmp(TraceList.dataSetName,'kikNet') );
idx_nga  = find(strcmp(TraceList.dataSetName,'ngawest1'));

scsnList = TraceList.selectSubList(idx_scsn);
jpList   = TraceList.selectSubList(idx_jp);
ngaList  = TraceList.selectSubList(idx_nga);

M  = [8,7,6,5,4,3,2];
nm = numel(M);
datTable = zeros(nm-1,6);

for im = 1:nm-1
   mup = M(im);
   mlo = M(im+1);
   
   datTable(im,1) = numel(find(jpList.m>=mlo & jpList.m<mup & jpList.hypDist< 25));
   datTable(im,2) = numel(find(jpList.m>=mlo & jpList.m<mup & jpList.hypDist>=25));
   
   datTable(im,3) = numel(find(ngaList.m>=mlo & ngaList.m<mup & ngaList.hypDist< 25));
   datTable(im,4) = numel(find(ngaList.m>=mlo & ngaList.m<mup & ngaList.hypDist>=25));

   datTable(im,5) = numel(find(scsnList.m>=mlo & scsnList.m<mup & scsnList.hypDist< 25));
   datTable(im,6) = numel(find(scsnList.m>=mlo & scsnList.m<mup & scsnList.hypDist>=25));

end

printmat(datTable, 'No. of records', sprintf('M%i ',M), 'JP_NS JP_FS NGA_NS NGA_FS SCSN_NS SCSN_FS' )

filename = 'dataTable.xlsx';
writetable(datTable,filename,'Sheet',1,'Range','D1')