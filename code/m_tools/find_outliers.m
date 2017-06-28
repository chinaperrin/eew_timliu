fprintf(1,['\nExtracting max-amplitudes in ',num2str(nr_snip*0.5),'sec window from ...\n'])

A = zeros(ntr_raw,nbands);
for iband = 1:nbands
    
    fprintf(1,['  ',num2str(iband),'th frequency band\n'])
    A(:,iband) = cellfun(@(x) x(iband,nr_snip), GlobalList.amax);
end

A_log = log10(A);
fprintf(1,'                                   DONE\n')


lgc_m   = GlobalList.m<=4;
lgc_a   = A_log(:,5)>=0;
lgc_tot = (lgc_m & lgc_a);

idx_outlier = find(lgc_tot);



nout = numel(idx_outlier);

for iout = 1:nout
    traceFullName = GlobalList.fullName{idx_outlier(iout)};
    [S,meta] = read_any_trace(traceFullName,0);
    
    gcf; clf; hold on
    plot(meta.t,S.raw)
    1+1;

end


% write ascii list of outliers
mycell        = cell(nout,2);

% filenames
splitList = regexp({GlobalList.fullName{idx_outlier}}','/','split');
splitList = [splitList{:}];
nsp       = numel(splitList);
tens      = linspace(10,nsp,nsp/10);
fileList  = splitList(tens)';

mycell(:,1) = fileList;
mycell(:,2) = {GlobalList.eqDate{idx_outlier}}';

[nrows,ncols] = size(mycell);

filename = 'outlier_list.txt';
fid = fopen(filename, 'w');

for row=1:nrows
    fprintf(fid, '%s %s \n', mycell{row,:});
end

fclose(fid);


