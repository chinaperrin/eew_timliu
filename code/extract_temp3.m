load out_tr137595_pid14439_all3275.mat;
zList = out.zList;
idx = find(cellfun(@(x) ~isempty(x), zList.var.v1));
procList = zList.selectSubList(idx);
[m, n] = size(procList.fullName);
currentfolder = pwd;
name = Q;  %replace Q with variable name
LEN = J;   %replace 'J' with number of variable subfields
cd '../ex_data'
mkdir(name)
cd(name)
idxFirstEst = cellfun(@(x) find(x.gba_ss.rh.hat2,1,'first'), procList.var.v1,'uniformOutput',0);
idxFirstEst = cell2mat(idxFirstEst);
for i = 0:LEN
    if i < 10
        filename = char(strcat(name, '_0', int2str(i), '.txt'));
    end
    if i > 9
        filename = char(strcat(name, '_', int2str(i), '.txt'));
    end
    fileID = fopen(filename, 'w');
    for x = 1:m
        if (idxFirstEst(x) + i) <= 40
            dat = K(idxFirstEst(x)+i);  %replace K with object sub-field
            if isfloat(dat)
                fprintf(fileID, '%.4f\n', dat);
            end
            if ischar(dat)
                fprintf(fileID, '%s\n', dat);
            end
        else
            fprintf(fileID, '0\n');
        end
    end
    fclose(fileID);
end
cd(currentfolder);