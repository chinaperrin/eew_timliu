function [outlist] = read_nga_metafiles_string(fName)

% Imports manually produced single columns from the NGA West1 Flatfile
% For details on these files, see /scratch/memeier/VS/data/wform/ngawest/meta/readme
% menandrin@gmail.com, 130912


fid    = fopen(fName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

eolLim      = '\r?\n';              % End of line
LinePattern = ['[^\r\n]*', eolLim]; % Line
thelines    = regexp(StrRay,LinePattern,'match');
outlist     = thelines';
%outlist     = cellfun(@(x) strtrim(x), outlist,'uniformOutput',0);
%outlist     = cellfun(@(x) x(1), thelines');

%idx_empty          = find(cellfun(@(x) isempty(x), outlist));
idx_empty          = find(cellfun(@(x) isempty(strtrim(x)), outlist));
outlist(idx_empty) = {'nofile/nofile'};