function save_padded_nga_wforms(s,srcFullName,targetFullName)

ucmd = strcat(['head -n4 ',srcFullName,' > ',targetFullName]);
unix(ucmd);


ncol = 5;
ns   = numel(s);


nln  = fix(ns/ncol);
% nrem = rem(ns,ncol);      % Neglect out last entries
sout = reshape(s(1:ncol*nln),ncol,nln)';

fid = fopen(targetFullName, 'a');
for iln=1:nln
    fprintf(fid, '%f %f %f %f %f', sout(iln,:));
    fprintf(fid, '\n');
end
fclose(fid);
