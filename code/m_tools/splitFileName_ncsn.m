function fNameParts = splitFileName_ncsn(fileFullName); 

% Separate path from rest of filename
slashIdx = regexp(fileFullName,'\/');
if ~isempty(slashIdx); fNameParts.path = fileFullName(1:slashIdx(end)); 
                       fileName        = fileFullName(slashIdx(end)+1:end);
else                   fNameParts.path = [];
                       fileName        = fileFullName;
end

pointIdx = regexp(fileName,'\.');

fNameParts.traceName  = fileName;
fNameParts.firstField = fileName(1:pointIdx(1)-1);
fNameParts.network    = fileName(pointIdx(1)+1:pointIdx(2)-1);
fNameParts.station    = fileName(pointIdx(2)+1:pointIdx(3)-1);
fNameParts.channel    = fileName(pointIdx(3)+1:pointIdx(4)-1);
fNameParts.extension  = fileName(pointIdx(end)+1:end);
fNameParts.recordName = fileName(1:end-7); 

if numel(pointIdx)==5; fNameParts.extrafield = fileName(pointIdx(4)+1:pointIdx(5)-1);
else                   fNameParts.extrafield = [];
end

