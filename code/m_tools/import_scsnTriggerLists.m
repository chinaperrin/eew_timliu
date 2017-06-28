function outTable = import_scsnTriggerLists(ListFileName)

fprintf(1,sprintf('Importing triggerList: %s ... ',ListFileName))
fid    = fopen(ListFileName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

codd = 0;

eolLim      = '\r?\n';                  % End of line
LinePattern = ['[^\r\n]*', eolLim];     % Line
thelines    = regexp(StrRay,LinePattern,'match');
nlines      = numel(thelines);

if (nlines==0); error('list could not be read.'); end

list.nw             = cell (nlines,1);
list.sta            = cell (nlines,1);
list.chan           = cell (nlines,1);
list.dateTimeString = cell (nlines,1);
list.subsec         = zeros(nlines,1);
list.amp1           = zeros(nlines,1);
list.amp2           = zeros(nlines,1);
list.amp3           = zeros(nlines,1);
list.comment1       = cell (nlines,1);
list.comment2       = cell (nlines,1);

skipThisLine = false(nlines,1);

for iline=1:nlines

    ctLine  = thelines{iline};
    fields  = regexp(ctLine,'\s+', 'split');
    
    if ~strcmp(fields{1},'Warning')
        
        list.nw{iline}   = fields{1};
        list.sta{iline}  = fields{2};
        list.chan{iline} = fields{4};
        
        dateTimeString             = strrep(fields{5}     ,'T',',');
        dateTimeString             = strrep(dateTimeString,'Z','' );
        ptIdx                      = regexp(dateTimeString,'\.');
        
        list.dateTimeString{iline} = dateTimeString(1:ptIdx-1);
        list.subsec(iline)         = str2double(dateTimeString(ptIdx:end));
        
        list.amp1(iline) = str2double(fields{6});
        list.amp2(iline) = str2double(fields{7});
        list.amp3(iline) = str2double(fields{8});
        
        % Read quoted comments
        idxQuote = regexp(ctLine,'''');
        nQuotes  = numel(idxQuote);
        if nQuotes==2;     list.comment1{iline} = ctLine(idxQuote(1)+1:idxQuote(2)-1);
        elseif nQuotes==4; list.comment1{iline} = ctLine(idxQuote(1)+1:idxQuote(2)-1);
                           list.comment2{iline} = ctLine(idxQuote(3)+1:idxQuote(4)-1);
        elseif nQuotes==0; list.comment1{iline} = 'no comments';
        else               fprintf(1,'Odd number of quotes found\n.')
            codd = codd+1;
        end
        
    else
        skipThisLine(iline) = true;
    end
end

if codd~=0; fprintf(1,sprintf('%i entries with odd quote number in list.\n',codd)); end

% Before Matlab 2015:
fdNames = fieldnames(list);
for ifield = 1:numel(fdNames)
  evalString = sprintf('outTable.%s = list.%s(~skipThisLine,:)',fdNames{ifield},fdNames{ifield});
  eval(evalString)
end
% Matlab 2015+: outTable = outTable(~skipThisLine,:);

fprintf(1,'done.\n')