function strng = string2fileName(strng)

% Remove the following characters from string
rmStrng = {' ',';'};
for i = 1:numel(rmStrng)
    strng = strrep(strng,rmStrng{i},'');
end

% Replace the following characters with underscore '_'
uscStrng = {':',','};
for i = 1:numel(uscStrng)
    strng = strrep(strng,uscStrng{i},'_');
end