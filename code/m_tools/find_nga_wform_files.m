function [rawFullName] = find_nga_wform_files(fname,dirname,wformDir,o_verbose)

% Searches for the file <fname> in <wformDir>, returns empty if not found
if (o_verbose)
    fprintf(1,['   Searching for ',fname,' by the name of ',fname,'* in dir ',dirname,'\n\t']);
end

[~,result]  = unix(['find ',wformDir, ' -name ', fname,'*']);
eolLim      = '\r?\n';              % End of line
LinePattern = ['[^\r\n]*', eolLim]; % Line
resultlines = regexp(result,LinePattern,'match');

% Find the one with the right directory name
match    = regexp(resultlines,dirname);
idxMatch = find(cellfun(@(x) ~isempty(x), match));

if (numel(idxMatch)==1)
    rawFullName = strtrim(resultlines{idxMatch});
    if (o_verbose); fprintf(1,['One match: ',rawFullName,'\n']); end
   
elseif (numel(idxMatch)>1)
    fprintf(1,['Multiple matches: ',resultlines{idxMatch},' \n\n\t WHATNOW?'])
    pause
    
elseif (numel(idxMatch)==0)
    if (o_verbose); fprintf(1,['No matches: ',result,'\n']); end
    rawFullName = [];
end