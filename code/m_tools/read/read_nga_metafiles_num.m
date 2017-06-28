function [outmat] = read_nga_metafiles_num(fName)

% Imports manually produced single columns from the NGA West1 Flatfile
% For details on these files, see /scratch/memeier/VS/data/wform/ngawest/meta/readme
% menandrin@gmail.com, 130912


fid    = fopen(fName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

eolLim      = '\r?\n';              % End of line
LinePattern = ['[^\r\n]*', eolLim]; % Line
thelines    = regexp(StrRay,LinePattern,'match');
outvect     = cellfun(@(x) str2num(x), thelines,'UniformOutput',0)';

% Write -1 into empty cells
idxEmpty          = find(cellfun(@(x) isempty(x),outvect));
outvect(idxEmpty) = {-1};
outmat            = cell2mat(outvect);