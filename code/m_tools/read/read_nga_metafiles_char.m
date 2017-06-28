function [outlist] = read_nga_metafiles_char(fName)

% Imports manually produced single columns from the NGA West1 Flatfile,
% where each column contains a single character.
% For details on these files, see /scratch/memeier/VS/data/wform/ngawest/meta/readme
% menandrin@gmail.com, 130912


fid    = fopen(fName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

eolLim      = '\r?\n';              % End of line
LinePattern = ['[^\r\n]*', eolLim]; % Line
thelines    = regexp(StrRay,LinePattern,'match');
outlist     = cellfun(@(x) x(1), thelines');