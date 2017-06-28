function cell2txt(myCell,fileName,)
[nrows,ncols] = size(myCell);

fid = fopen(fileName, 'w');

for row=1:nrows
    fprintf(fid, '%s\n', mycell{row,:});
end 
 fclose(fid);
