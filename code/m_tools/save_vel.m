function save_vel(s,traceFileName,traceName,method)


% Outfile name
slashPositions = find(traceFileName == '/');
velDirName     = strcat(traceFileName(1:slashPositions(end)),'vel/');
if (~exist(velDirName,'dir'))
    ucmd = strcat(['mkdir ',velDirName]); unix(ucmd);
end
pointPositions = find(traceName == '.');
velFileName = strcat(velDirName,traceName(1:pointPositions(end)-1),'.VEL');


if (method == 1)
    tic
    dlmwrite(velFileName,s)
    toc
    
elseif (method == 2)
    velFileName2 = strcat(velDirName,traceName(1:pointPositions(end)-1),'.VEL2');
    tic;
    fid = fopen(velFileName2, 'w+');
    for i=1:size(s,1)
        fprintf(fid, '%f ', s(i));
        fprintf(fid, '\n');
    end
    fclose(fid);
    toc
    
elseif (method == 3)
    velFileName3 = strcat(velDirName,traceName(1:pointPositions(end)-1),'.VEL3');
    tic;
    csvwrite(velFileName3, s);
    toc;
end