function [status] = scp_wform(traceFullName,scpServerName)

if nargin<2; scpServerName='trinitite'; end

slashIdx = regexp(traceFullName,'/');
pathName = traceFullName(1:slashIdx(end));
 
if (~exist(pathName,'dir')); 
    fprintf(1,['Creating directory\n  ',pathName,'\n'])
    ucmd = strcat(['mkdir -p ',pathName]); 
    unix(ucmd); 
end

ucmd = sprintf('scp -rp memeier@%s.gps.caltech.edu:%s %s',scpServerName,traceFullName,pathName);
%ucmd = strcat(['scp -rp memeier@trinitite.gps.caltech.edu:',traceFullName, ' ',pathName]);
%ucmd = strcat(['scp -rp men@bigstar02.ethz.ch:',traceFullName, ' ',pathName]);
[status,~] = unix(ucmd);