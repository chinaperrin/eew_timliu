function [fullNameZ,fullNameE,fullNameN] = find_3corecs_fullNames(traceFullName)
% Find and read corecorded traces of trace called <traceFullName>. Looks
% for corecs only in same directory as <traceFullName>. Does not need/use
% an existing traceList. Only works for SCSN sac-files.
%
% menandrin@gmail.com, 140309


% Find out what sort of waveform it is ------------------------------------
slashIdx       = regexp(traceFullName,'/');
pathName       = traceFullName(1:slashIdx(end));
ptIdx          = regexp(traceFullName,'\.');
extension      = traceFullName(ptIdx(end):end);
jp_extensions  = {'.EW','.EW1','.EW2','.NS','.NS1','.NS2','.UD','.UD1','.UD2'};


if strcmp(extension,'.sac')       % - - - - - - - - -   SOUTHERN CALIFORNIA

    recordFullName = traceFullName(1:end-5);
    fullNameZ      = strcat(recordFullName,'Z.sac');
    fullNameE      = strcat(recordFullName,'E.sac');
    fullNameN      = strcat(recordFullName,'N.sac');
    
    
elseif ismember(extension,jp_extensions)    % - - - - - - - - - -     JAPAN
    
    orntString = traceFullName(ptIdx(end):ptIdx(end)+2);
    fullNameZ  = strrep(traceFullName,orntString,'.UD');
    fullNameE  = strrep(traceFullName,orntString,'.EW');
    fullNameN  = strrep(traceFullName,orntString,'.NS');
    
    
elseif strcmp(extension,'.AT2')             % - - - - - - - - -   NGA West1
    
    uscIdx         = regexp(traceFullName,'_');
    recordFullName = traceFullName(1:uscIdx(end));
    fileList       = dir([recordFullName,'*']);
    nrec           = size(fileList,1);
    
    if nrec==3
        
        % Find vertical one
        names   = {fileList(1).name; fileList(2).name; fileList(3).name};
        idxVert = find(cellfun(@(x) ~isempty(x),regexp(names,'_V.AT2')));
        idxHor  = setdiff(1:3,idxVert);
        
        fullNameZ = strcat([pathName,names{idxVert  }]);   % Assign vertical one
        fullNameE = strcat([pathName,names{idxHor(1)}]);   % Randomly assign horizntal ones
        fullNameN = strcat([pathName,names{idxHor(2)}]);
        fprintf(1,'8ung: horizontal components are randomly assigned.\n')

    else
        
        fprintf(1,sprintf('only %i instead of 3 components have been found',nrec))
        fullNameZ = [];
        fullNameE = [];
        fullNameN = [];
    end
else
    fprintf(1,'EXTENSION NOT FOUND, WHAT NOW?\n')
    pause
end