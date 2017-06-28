function cp_wform_lacie(traceFullName,pathName)
% Copy the waveform from external HD to laptop HD

if (~exist(pathName,'dir')); 
    fprintf(1,['Creating directory\n  ',pathName,'\n'])
    ucmd = strcat(['mkdir ',pathName]); 
    unix(ucmd); 
end

ucmd = strcat(['cp -rp /Volumes/LaCie',traceFullName, ' ',pathName]);
unix(ucmd);