load_traceList_basic;
[m, n] = size(zList.fullName);
currentfolder = pwd;
name = Q;  %replace Q with variable name
LEN = J;   %replace 'J' with number of variable subfields
cd '../../ex_data'
mkdir(name)
cd(name)
for i = 1:LEN
    if i < 10
        filename = char(strcat(name, '_0', int2str(i), '.txt'));
    end
    if i > 9
        filename = char(strcat(name, '_', int2str(i), '.txt'));
    end
    fileID = fopen(filename, 'w');
    for x = 1:m
        dat = K(i);  %replace K with object sub-field
        if isfloat(dat)
            fprintf(fileID, '%.4f\n', dat);
        end
        if ischar(dat)
            fprintf(fileID, '%s\n', dat);
        end
    end
    fclose(fileID);
end
cd(currentfolder);