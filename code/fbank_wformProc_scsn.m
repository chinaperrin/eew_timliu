function fbank_wformProc_scsn(wformDir); dbstop if error
%fbank_wformProc_scsn('/scratch/memeier/data/socal/M6p/');

% Goes through specified stp-waveform directory, through all event dubdirs, 
% processes all sacfiles of event subdirs which have a *.phase file. 
% 1) adds all traces with hypocentral distances < hD_max to a traceList 
%    object, and all other traces to skipList object
% 2) goes through all entries in traceList, tests them for several quality 
%    requirements (e.g. minimum sampling rates); if all requirements are 
%    fulfilled, trace is processed through filter bank; if not trace is 
%    added to skipList.
%
% menandrin@gmail.com, 131104

%clear all; wformDir  = '/scratch/memeier/data/socal/M6p/';


addpath(genpath('../'))
addpath(genpath('../../'))
addpath(genpath('../../../../matlab/'))

global fOrder fMode ftSize

configFileName = 'cfg/config_wformProc_scsn_GbA_i39.m';

% Mark NoCal events such that correct units are used (NoCal use SI units,
% SoCal use cm)
noCalList = {};
%noCalList = {'72282711';'999'};
if ~isempty(noCalList)
    fprintf(1,'Events in NoCal-List:\n')
    fprintf(1,sprintf('\t.%s\n',noCalList{:}))
    fprintf(1,'Are all NoCal-events included in List?')
    pause
end

%wformDir   = '/scratch/memeier/finder/wforms/cali/';
%wformDir   = '/scratch/memeier/data/shakealert/mam/';
%wformDir   = '/scratch/memeier/data/shakealert/tests/mam_test_CVS/';
%wformDir   = '/scratch/memeier/data/shakealert/tests/mam_test3/';
%wformDir  = '/scratch/memeier/data/socal/M6p/';
%wformDir  = '/scratch/memeier/data/socal/M5/';
%wformDir  = '/scratch/memeier/data/socal/M4/';
%wformDir  = '/scratch/memeier/data/socal/testM4/';
%wformDir  = '/scratch/memeier/data/socal/borregoSprings16/';
%wformDir  = '/scratch/memeier/data/socal/elizabeth/';
%dbstop if error

 


%% Overview
%
%  0. Parameter definition
%
%  1. Load traceList or create one by going through all event directories.
%     Select traces based on criteria that can be determined without
%     actually reading the waveform, e.g. hypocentral distance, existence
%     of phase picks in *.phase files, etc.
%
%  2. Go through traceList, for each entry ...
%
%       2.1 Read trace-info 
%           Quality control (QC) 1: maxCounts, no. of corecs
%           Find vertical corec
%
%       2.2 Load waveform
%           QC 2: max. & min. signal length, minimum SR
%
%       2.3 p-Pick (index)
%           Load vertical component
%           If scsn-network pick available, use it, else run picker on vertical wform
%           Measure SNR between before & after pick
%           QC 3: too high or too low SNR?
%
%       2.4 Evalute pick time (as opposed to index)
%           QC 4: sensible pick time?
%
%       2.5 Check for outliers
%           Highpass filter trace, with <fLow_prefilt>
%           Differentiate/integrate to acc/vel/dsp
%           Measure PGA/PGV/PGD 
%           Measure nbpGA/nbpGV/nbpGD in all subbands
%           Check for approximate consistency of PGV with CH2007
%
%       2.6 Pass waveform through filterBank 
%           ==> maximum velocities = fct(freq,time) <==
%           Measure & save pre-pick noise on velocity trace on full wform 
%           Plot results
%           Move all skipped traces to skipList
%           Save output
 

%% Notes
%  - TIME: for each time, read_sac_trace3.m generates a t-vector <t> with
%    t(t0) = 0. For each trace, store the first value <ts> of this time 
%    vector, the sampling rate, as well as the absolute origin time/date. 
%  - out of pure laziness for stations from longList, I assume stAltitude = 0m
%  - the saved snr is always on the velocity traces, but what i use for
%    keeping or discarding traces (snr>!20) is on acc trace for SM records
%  - when reading the SCSN sac files I have to round the sampling interval
%    to get integer sampling rates.



%% PARAMETERS AND SETTINGS

% READ CONFIG FILE ........................................................
fileId = fopen(configFileName);
StrRay = fscanf(fileId,'%c');       % Read the entire text file into
fclose(fileId);                     % an array

eolLim      = '\r?\n';                   % End of line
LinePattern = ['[^\r\n]*', eolLim];      % Line
thelines    = regexp(StrRay,LinePattern,'match');
nl          = numel(thelines);
for il = 1:nl
    oneLine = thelines{il};
    eval(oneLine)
end


% READ EXTERNAL FILES .....................................................
load(blackListFullName)
load(pxListFullName)
load(finSrcFullName);

stationList_cali     = import_stationlist_fB(stFileName);
stationLongList_cali = import_stationLongList(stFileName2);


% SETTING PATH- FILE AND DIRECTORY NAMES ..................................
wformDir1   = strcat(wformDir,'cms/');
dirAppendix = strcat(['_',strrep(num2str(snpLength),'.','p')]);
if     strcmp(fMode,'causal' ); outDir = strcat(['out/i',num2str(iN),'/causal' ,dirAppendix,dirAppendix2,'/']);
elseif strcmp(fMode,'acausal'); outDir = strcat(['out/i',num2str(iN),'/acausal',dirAppendix,dirAppendix2,'/']);
end

TraceListFullName = strcat([wformDir,outDir,'trList.mat']);
SkipListFullName  = strcat([wformDir,outDir,'skipList.mat']);
outFileName       = strcat([wformDir,outDir,'fbOut.mat']);
logFileName       = strcat([wformDir,outDir,'fB_log.txt']);
diaryFileName     = strcat([wformDir,outDir,'diary.txt']);

if (~exist(strcat([wformDir,outDir]),'dir'))
    unix(['mkdir -p ',strcat([wformDir,outDir])]);
else
    if (exist(outFileName,'file') & (o.saveOut) );
        fprintf(1,'WARNING: outfile already exists. If you proceed, all existing file will be overwritten! Think twice.')
        pause
    end
end


if o.saveOut
    currentDir    = pwd;
    [~,starttime] = unix('date');
    scrpt         = mfilename('fullpath'); scrpt = strcat([scrpt, '.m']);
    unix(['echo "This traceList was compiled with the following script on ',starttime,'" >  ' logFileName]);
    unix(['echo "  executed in ',currentDir,'" >>  ' logFileName]);
    unix(['cat ', scrpt, ' >> ', logFileName]);
    
    if (exist(diaryFileName,'file')); ucmd = strcat(['rm ',diaryFileName]); unix(ucmd); end
    diary(diaryFileName)
end










%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%      A. Load or gather list of all traces that have a broadband P-pick  %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Initialise global traceList
[~,ntrMax] = unix(['find ', wformDir1,' -name *.sac | wc -l']); 
%[~,ntrMax] = unix(['find ', wformDir1,' -name *',channelPattern,'.sac | wc -l']);
ntrMax     = str2num(ntrMax);

TraceList  = traceList(ntrMax);
traceIdx   = 0;
flg_tooFar = false(ntrMax,1);

fprintf(1,'\n\n\n\n-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n')
fprintf(1,['Compiling a list of all traces for the specified data set ... \n', ...
    '   --> ',num2str(ntrMax),' matching sac-files found ... \n'])

% Get a list of sub-directories (=events) in the specified wform-directory
dirList = dir(wformDir1);
eqIdx   = cellfun(@(x) str_isnumeric(x), {dirList.name}');
eqList  = dirList(logical(eqIdx));
neq     = size(eqList,1);
fprintf(1,[' ',num2str(neq),' waveform directories / events:\n\n'])

% For each sub-directory (= event) ...
for ieq = 1:neq
    
    hasFinSrcModel = false;
    
    eqName        = getfield(eqList,{ieq}, 'name');
    fileList      = dir(strcat(wformDir1,eqName));                % get a list of all files
    fileList      = {fileList.name}';
    phaseFileName = strcat(wformDir1,eqName,'/',eqName,'.phase'); % find pickfile
    
    fprintf(1,['\n ... event directory ',num2str(ieq),' out of ',num2str(neq),': ',eqName,' ... '])
    
    if ~ismember(eqName,noCalList); dsName = 'scsn';
    else                            dsName = 'ncsn';
    end
    
    % If there is a phaseFile ...
    if (exist(phaseFileName,'file'))
        
        % Find all sac-files with specified channel description
        sacIdx       = regexp(fileList,'.sac');
        isSac        = cellfun(@(x) ~isempty(x), sacIdx);
        chanIdx      = regexp(fileList,channelPattern);
        rightChannel = cellfun(@(x) ~isempty(x), chanIdx);
        sacList      = fileList(isSac &rightChannel);
        nsac         = numel(sacList);
        %sacIdx = regexp(fileList,strcat(channelPattern,'.sac')); --> does not work for NCSN files with station code or double dots
        
        % Read phase-file
        % Units:    eqZ in [km]; stAlt in [m]
        [eqType,eqDate,t0,eqLat,eqLon,eqZ,magn,magnType,phasePxList] = read_phase_file(phaseFileName);
        %[eqType,eqDate,t0,eqLat,eqLon,eqZ,magn,magnType,phasePxList] = read_phase_file_prei37(phaseFileName);
        
        % Write into eventList
        %eventIdx = eventIdx + 1;
        %evtList.addEvent(eventIdx,9999,eqLat,eqLon,eqZ,eqDate,magn,0);
        % Find finite source model
        if magn>=6
            ifault = find(strcmp(eqName,allFaults.eqname));
            if ~isempty(ifault) &~isempty(allFaults.segments{ifault}.lat)
                thisFault      = allFaults.segments{ifault};
                npoints        = numel(thisFault.lat);
                hasFinSrcModel = true;
            end
        end
        
        % Find maximum distance for which records should be considered
        if o.mdeprmax; maxFltDist = get_max_dist_with_gmPrime(magn,pgaPrime);
                       if maxFltDist<200; maxFltDist=200; end
        else           maxFltDist = maxFltDistFix;
        end

        ct1=0;
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Process all sac-traces in event-directory
        for isac = 1:nsac
            
            traceIdx      = traceIdx + 1;
            traceName     = sacList{isac};
            traceFullName = strcat(wformDir1,eqName,'/',traceName);
            
            nmParts       = splitFileName_ncsn(traceFullName);
            stationName   = nmParts.station;
            networkName   = nmParts.network;
            % Find vertical p-pick for this record, preferentially from a broadband channel
            %pointIdx        = regexp(traceName,'\.');
            %stationName     = traceName(pointIdx(2)+1:pointIdx(3)-1);
            %networkName     = traceName(pointIdx(1)+1:pointIdx(2)-1);
            [has_pPx,allPx] = get_pPick(phasePxList,networkName,stationName);
            
            % ... if there is one, add all corecorded traces to traceList
            if (has_pPx)
                TraceList.px.p.hasPx(traceIdx)  = true;
                TraceList.dataSetName{traceIdx} = sprintf('%sPx',dsName);
                %TraceList.dataSetName{traceIdx} = 'scsnPx';
                stLat = allPx.stLat;
                stLon = allPx.stLon;
                stAlt = allPx.stAlt;    % [m]
                tppx  = allPx.tppx;
            else
                
                TraceList.px.p.hasPx(traceIdx) = false;
                TraceList.dataSetName{traceIdx} = dsName;
                %TraceList.dataSetName{traceIdx} = 'scsn';
                
                [stLat,stLon,networkName,stationName,chan] = get_stationCoords_cali(traceFullName,stationLongList_cali);
                
                stAlt = 0;
                tppx  = 0;
                
                % If not found, set stCoords such that hypDist will be too large and trace will be rejected
                if isempty(stLat); stLat = 0;
                    stLon = 0;
                end
            end
            
            % Compute source/station distances
            hypDist = hypoDistance(eqLat,eqLon,eqZ*1e3,stLat,stLon,-stAlt);   % Input in [m] & [deg], output in [m]
            epiDist = hypoDistance(eqLat,eqLon,0,      stLat,stLon,-stAlt);
            hypDist = hypDist/1e3;                                            % Save as km
            epiDist = epiDist/1e3;                                            % Save as km
            
            if hasFinSrcModel
                fdVect = zeros(npoints,1);
                for ipt=1:npoints
                    fdist_m     = hypoDistance(thisFault.lat(ipt),thisFault.lon(ipt),thisFault.z(ipt)*1e3,stLat,stLon,-stAlt);   % Input in [m] & [deg], output in [m]
                    fdVect(ipt) = fdist_m*1e-3;
                end
                fltDist = min(fdVect);
            else
                fltDist = hypDist;
            end
            
            if fltDist>maxFltDist; flg_tooFar(traceIdx)=true; ct1=ct1+1; end

            
            % Save all relevant info
            TraceList.fullName{traceIdx}      = traceFullName;
            TraceList.eq.t0{traceIdx}         = t0;
            TraceList.eq.date{traceIdx}       = eqDate;
            TraceList.eq.name{traceIdx}       = eqName;
            TraceList.eq.m(traceIdx)          = magn;
            TraceList.eq.mType{traceIdx}      = magnType;
            TraceList.eq.lat(traceIdx)        = eqLat;
            TraceList.eq.lon(traceIdx)        = eqLon;
            TraceList.eq.z(traceIdx)          = eqZ;
            
            TraceList.dist.hyp(traceIdx)      = hypDist;
            TraceList.dist.epi(traceIdx)      = epiDist;
            TraceList.dist.flt(traceIdx)      = fltDist;
            
            TraceList.px.p.t(traceIdx)        = tppx;         % Always use broadband P-pick
            
            TraceList.station.name{traceIdx}  = stationName;
            TraceList.station.nw{traceIdx}    = networkName;
            TraceList.station.lat(traceIdx)   = stLat;
            TraceList.station.lon(traceIdx)   = stLon;
            TraceList.station.alt(traceIdx)   = stAlt;
            
            TraceList.station.bcode{traceIdx} = nmParts.channel(1);
            TraceList.station.icode{traceIdx} = nmParts.channel(2);
            TraceList.station.ocode{traceIdx} = nmParts.channel(3);
            
            %             fprintf(1,sprintf('     %i/%i traces, M%3.1f at R(e/h/f) = %5.1f/%5.1f/%5.1fkm',isac,nsac,magn,epiDist,hypDist,fltDist))
            %             if flg_tooFar(traceIdx); fprintf(1,'\t too far.\n')
            %             else                     fprintf(1,'\n')
            %             end

        end
        fprintf(1,'%i/%i records were too far.\n',ct1,nsac)
        
    else
        fprintf(1,'No phase-file found, skipping all traces of this event\n\n')
    end
end

%if o.saveOut; save('scsn_after0stLoop.mat'); end
%save('scsn_after0stLoop.mat')

% Throw out those with too large hypocentral distances
%flg_tooFar = ( (TraceList.dist.hyp>hD_max &TraceList.eq.m<Mthd) |(TraceList.dist.hyp>hD_max2 &TraceList.eq.m>=Mthd) );
%flg_tooFar                    = (TraceList.dist.hyp>hD_max);
idx_tooFar                    = find( flg_tooFar);
idx_closeEnough               = find(~flg_tooFar);
skipReason                    = 'too far';
TraceList.comment(flg_tooFar) = {skipReason};                           % XXXXXXX not sure how this one is to be done
fprintf(1,'check me.\n')

% Add all traces with flg_tooFar(itr) = true to SkipList
fprintf(1,['    adding ',num2str(numel(idx_tooFar)),' trace to SkipList ...'])
SkipList  = TraceList.selectSubList(idx_tooFar);
TraceList = TraceList.selectSubList(idx_closeEnough);
fprintf(1,['    done.\n'])  % Remove skipped traces and outliers from traceList

% If there are wform-directories that do not contain a *.phase-file, the above
% loop will skip them, which leads to empty list-entries at the end of
% the TraceList. Get rid of the empty lines:
%TraceList.consolidate
idxNotEmpty = find(TraceList.eq.m~=0);
TraceList   = TraceList.selectSubList(idxNotEmpty);



TraceList.prop.o                 = o;
TraceList.prop.fc                = fc;
TraceList.prop.snpLength         = snpLength;
TraceList.prop.outDirName        = strcat([wformDir,outDir]);
TraceList.prop.envShift          = envShift;
TraceList.prop.envWidth          = envWidth;
TraceList.prop.iN                = iN;
TraceList.prop.px                = px;
TraceList.prop.blackListFullName = blackListFullName;
TraceList.prop.pxListFullName    = pxListFullName;
TraceList.prop.prePxWindow       = prePxWindow;
TraceList.prop.tnoise            = tnoise;
TraceList.prop.maxCounts         = maxCounts;
TraceList.prop.te_onset          = te_onset;
TraceList.prop.tpsgap            = tpsgap;
TraceList.prop.accThresholds     = accThresholds;
TraceList.prop.velThresholds     = velThresholds;
TraceList.prop.immiThresholds    = immiThresholds;

TraceList.prop.azi.bazIntervalVector = bazIntervalVector;

TraceList.prop.gmThresh.mmi.tex  = immiThresholds;
TraceList.prop.gmThresh.acc.tex  = accThresholds;
TraceList.prop.gmThresh.vel.tex  = velThresholds;

TraceList.prop.rmFieldList             = rmFieldList;
TraceList.prop.files.wformDir          = wformDir;
TraceList.prop.files.outDir            = outDir;
TraceList.prop.files.logFileName       = logFileName;
TraceList.prop.files.diaryFileName     = diaryFileName;
TraceList.prop.files.TraceListFullName = TraceListFullName;
TraceList.prop.files.SkipListFullName  = SkipListFullName;
TraceList.prop.files.outFileName       = outFileName;
TraceList.prop.files.configFileName    = configFileName;
TraceList.prop.files.configLines       = thelines;

if o.saveOut
    fprintf(1,'\n     Saving TraceList ... ')
    save(TraceListFullName,'TraceList');
    save(SkipListFullName,'SkipList');
    fprintf(1,' done\n')
end


% Print trace summary
ntr  = size(TraceList.eq.m,1);
nskp = size(SkipList.eq.m,1);
summary_trList(TraceList,SkipList)


%                '-._   ```"""---.._
%             ,-----.:___           `\  ,;;;,
%              '-.._     ```"""--.._  |,%%%%%%              _
%              ,    '.              `\;;;;  -\      _    _.'/\  
%            .' `-.__ \            ,;;;;" .__{=====/_)==:_  ||  - - -  support the filter bank initiative!
%       ,===/        ```";,,,,,,,;;;;;'`-./.____,'/ /     '.\/
%      '---/              ';;;;;;;;'      `--.._.' /
%     ,===/                          '-.        `\/
%    '---/                            ,'`.        |
%       ;                        __.-'    \     ,'
%       \______,,.....------'''``          `---`
% 
if (o.silly); print_asciiart('m_tools/asciiart/angel.txt'); end %pause(.5)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%       B0. Sort TraceList such that corecs are next to each other        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% With station(?)-codes in filenames (__E.01.sac) different components of
% the same station can be located at some distance from each other. This
% trips up the next section (B.) where corecs are expected to be located
% next to each other in TraceList. Make sure that's the case.
% Also: make sure vertical comes in third place such that tOcc will be
% saved with vertical
if o.sortTraceList

    % Assign traceList entries to earthquakes
    eqs = compile_eqs_structure(TraceList,[]);
    neq = numel(eqs.m);
    
    % For each earthquake, make sure co-recorded traces are next to each other
    %fprintf(1,'Making sure that corecorded traces are next to each other in traceList ... ')
    sortedList = traceList(0);
    
    for ieq=1:neq
   
        fprintf(1,sprintf('Sorting %i/%ieqs ... ',ieq,neq))
        %print_iteration_numbers(ieq,neq,'tens')
            
        thisList = TraceList.selectSubList(eqs.traceId{ieq});
        ntr      = numel(thisList.eq.m);
        itr      = 1;
        
        sortIdx     = [];
        dropIdx     = [];
        hasBeenProc = false(ntr,1);
        
        while itr<=ntr
            
            if ~hasBeenProc(itr);
                
                traceFullName = thisList.fullName{itr};
                traceSpecs    = textscan(traceFullName,'%s','delimiter','/');
                traceName     = traceSpecs{1}{end}; %recordName    = traceName(1:end-5);
                
                % Identify corecorded traces by searching for identical file-names with
                % different orientation-code
                chanIdx               = regexp(traceName,channelPattern); %traceName(chanIdx:chanIdx+2)
                orntIndex             = chanIdx+2;
                orntIndex             = orntIndex(end); % If station name also matches channelPattern ...
                searchName            = traceName;
                searchName(orntIndex) = '.';
                hits                  = regexp(thisList.fullName,searchName);
                idxCorecs             = find(cellfun(@(x) ~isempty(x), hits));
                nCorec                = numel(idxCorecs);
                thisList.var.v8(idxCorecs) = {searchName};
                
                if nCorec==3
                    zHit    = find(strcmp(thisList.station.ocode(idxCorecs),'Z'));
                    vertIdx = idxCorecs(zHit);
                    %sortIdx = [sortIdx; vertIdx; idxCorecs(idxCorecs~=vertIdx)];
                    sortIdx = [sortIdx; idxCorecs(idxCorecs~=vertIdx); vertIdx];
                else
                    dropIdx = [dropIdx; idxCorecs];
                end
                hasBeenProc(idxCorecs) = true;
            end
            itr=itr+1;
        end
        
        % Add to sortedList and skipList
        skipReason                 = 'not 3 corecs';
        thisList.comment(dropIdx) = {skipReason};
        ndrop = numel(dropIdx);
        if ndrop>0; fprintf(1,sprintf(' skipping %i/%i, no 3 corecs found.\n',ndrop,ntr))
        else        fprintf(1,'\n')
        end
        thisSkipList = thisList.selectSubList(dropIdx);
        SkipList.appendList(thisSkipList);

        tmpList    = thisList.selectSubList(sortIdx);
        sortedList.appendList(tmpList);
        
        clear tmpList thisSkipList
    end
    fprintf(1,'done.\n')
end
TraceList = sortedList;
ntr       = numel(TraceList.eq.m);
summary_trList(TraceList,SkipList)






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%          B. Process all traces of list                                  %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

fprintf(1,'\n\n\n-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -\n')
fprintf(1,'Running traces through filter bank ...\n\n')

flg_skipNow       = false(ntr,1);  % If this flag is raised, processing is aborted for this trace
flg_skipAfterProc = false(ntr,1);  % If this flag is raised, processing is continued and trace is added to SkipList at the end

hasScsnPx             = TraceList.px.p.hasPx;             % (the hasPpx will later be modified)
%TraceList.px.p.hasPx  = false(numel(TraceList.eq.m),1);   % This flag is used for indicating whehter a px
                                                          % has already been processed or not

% Make a backup copy of the uprocessed TraceList; the original list will be
% overwritten with the processed one.
% if (o.saveOut)
%     fprintf('Backin up unprocessed traceList\n')    
%     TraceListFullName_mod = strrep(TraceListFullName,'List.mat','List_unproc.mat');
%     unix(['cp ',TraceListFullName,' ',TraceListFullName_mod]);
% end

% Records come in "triplets", 3 components of each station. tripletIdx is 
% number of the record per triplet that currently being processed. 
% tripletIdx=0 is used as an indication that a new triplet should be started
% Note: if you comment out skipNow-flag, don't set tripletIdx to zero!
tripletIdx = 0;
ctref      = 0;
cvert      = 0;
fnameList  = [];

% For all traces in list ...
ntr_bak = ntr;
if ntr>ntrMaxAbs; ntr=ntrMaxAbs; end
%ntr = ntr-rem(ntr,3); % Make sure ntr is a multiple of three

for itr = 1:ntr
    
    clear t ns ppxIdx ppxIdx_z
    commentsSoFar = TraceList.comment{itr};
    if isempty(commentsSoFar); skipReasonList = {};
    else                       skipReasonList = {};
                               skipReasonList = [skipReasonList; commentsSoFar];
    end
    
    traceFullName = TraceList.fullName{itr};
    traceSpecs    = textscan(traceFullName,'%s','delimiter','/');
    traceName     = traceSpecs{1}{end};
    
    dsName = TraceList.dataSetName{itr};
    if regexp(dsName,'scsn'); isSocalTrace = true;
    else                      isSocalTrace = false;
    end
    
    % Check if record is black-listed
    searchName   = TraceList.var.v8{itr};
    idxBlackList = find(cellfun(@(x) ~isempty(x), regexp(blackList.fullName,searchName)));
    if ~isempty(idxBlackList) & ~flg_skipNow(itr)   % second condition is to avoid checking again with co-recorded traces
        hits                   = regexp(TraceList.fullName,searchName);
        idxCorecs              = find(cellfun(@(x) ~isempty(x), hits));
        flg_skipNow(idxCorecs) = true;
        thisSkipReason         = blackList.comment{idxBlackList};
        skipReasonList         = [skipReasonList; thisSkipReason];
    end

    if ~flg_skipNow(itr)
        
        % If itr is first of new triplet, identify corecorded traces
        if tripletIdx==0
            
            tripletIdx = 1;         % First trace of a triple
            
            % Clear out some values from previous trace triplet
            clear tOcc vertFullName outpx srawpx tz nsz srpx tzs

            % Identify location of corecorded traces in TraceLists
            searchName = TraceList.var.v8{itr};
            hits       = regexp(TraceList.fullName,searchName);
            idxCorecs  = find(cellfun(@(x) ~isempty(x), hits));
            zHit       = find(strcmp(TraceList.station.ocode(idxCorecs),'Z'));
            eHit       = find(strcmp(TraceList.station.ocode(idxCorecs),'E'));
            nHit       = find(strcmp(TraceList.station.ocode(idxCorecs),'N'));
            vertIdx    = idxCorecs(zHit);
            
            %             fNameParts = splitFileName_ncsn(traceFullName)
            %             if ~isempty(fNameParts.extrafield)
            %                 1+1;
            %             end
            
            % Quality control 1     % THESE CASES HAVE ALREADY BEEN FILTERED OUT IN SECTOIN B0.
            % -----------------
            %             if (nCorec<3)
            %                 flg_skipNow(idxCorecs) = 1;
            %                 thisSkipReason = 'less than 3 corecs';
            %                 skipReasonList = [skipReasonList; thisSkipReason];
            %                 %skipReason              = 'less than 3 corecs';
            %                 fprintf(1,[traceName, '\n\tskipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
            %                 %tripletIdx = 0;
            %             end
            %             if (nCorec>3)
            %                 flg_skipNow(idxCorecs) = 1;
            %                 thisSkipReason = 'more than 3 corecs';
            %                 skipReasonList = [skipReasonList; thisSkipReason];
            %                 %skipReason              = 'more than 3 corecs';
            %                 fprintf(1,[traceName, '\n\tskipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
            %                 %tripletIdx = 0;
            %             end
            %             if (isempty(vertIdx))
            %                 fprintf(1,'\n\t\t\thas already been skipped because there are no three traces ... ')
            %             end
            
            countsTraceFullName = strrep(traceFullName,'/cms/','/cnts/');
            hasCountsFile       = exist(countsTraceFullName,'file');
            if hasCountsFile==0
                flg_skipNow(idxCorecs) = 1;
                thisSkipReason = 'no counts file found for clipping check';
                skipReasonList = [skipReasonList; thisSkipReason];
                %skipReason              = 'no counts file found for clipping check';
                fprintf(1,[traceName, '\n\tskipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
                %tripletIdx = 0;
            end
        else
            tripletIdx = tripletIdx +1;
        end
    else
        tripletIdx=0;
    end

    fprintf(1,sprintf('\n%i/%i -- %i: %s ',itr,ntr,tripletIdx,traceName))
    

            
    % -----------------------------
    % B1. Read info from traceList 
    % -----------------------------
    if (~flg_skipNow(itr))
        
        % Load meta-info
        eqLat           = TraceList.eq.lat(itr);
        eqLon           = TraceList.eq.lon(itr);
        eqZ             = TraceList.eq.z(itr);
        eqYr            = TraceList.eq.date{itr};
        torg            = TraceList.eq.t0{itr};
        mtype           = TraceList.eq.mType{itr};
        hD              = TraceList.dist.hyp(itr);
        stLat           = TraceList.station.lat(itr);
        stLon           = TraceList.station.lon(itr);
        instrumentCode  = TraceList.station.icode{itr};
        orientationCode = TraceList.station.ocode{itr};
        
        isSM = (strcmp(instrumentCode,'L') ||strcmp(instrumentCode,'N') ||strcmp(instrumentCode,'G'));
        
        % Rough signal duration considerations: at 100km, s-arrival ~30s after 
        % t0; path duration ~10sec, source duration for m8.0 ~40s --> 80s
        m = TraceList.eq.m(itr);
        if m>=0; tmax = 20;  end
        if m>=5; tmax = 40;  end
        if m>=6; tmax = 60;  end
        if m>=7; tmax = 120; end

        
        % Quality control 2
        % -----------------        
        % Read uncorrected trace [counts] to check if it is clipped
        [out]                  = read_sac_trace3(countsTraceFullName);
        scounts                = out.sraw;
        tcounts                = out.t;
        maxObsCounts           = max(abs(scounts));
        %[scounts,~,~]         = read_sac_trace(countsTraceFullName);
        %[scounts,~,~,flgIssue] = read_sac_trace2(countsTraceFullName);
        
        %TMPTMPTMPTMPTMPTMP Compare phase file origin time with that of sac file
        sacdate          = strcat([TraceList.eq.date{itr},' ',torg]);
        tref_equals_torg = strcmp(sacdate,out.tref);
        if ~tref_equals_torg
            ctref = ctref+1;
            TraceList.comment(idxCorecs) = {'tref not equal t0'};
            fprintf(1,'\n\n\n\n\n\n\t\t  xxxxxxxxx    ttref not equal t0  xxxxxxx\n\n\n\n')
            %pause
        end    
        %TMPTMPTMPTMPTMPTMPTMPTMPTMPTMPTMPTMPTMPTMPTMPTMPTMPTMPTMPTMPTMPTMP
        
        % Skip traces ...
        if (maxObsCounts > maxCounts)                 % Clip test
            %flg_skipNow(idxCorecs) = 1;
            %flg_skipAfterProc(idxCorecs) = 1;
            thisSkipReason = 'clipped';
            skipReasonList = [skipReasonList; thisSkipReason];
            %TraceList.comment(idxCorecs) = {skipReason};
            fprintf(1,[traceName, '\n\tNOT skipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
        end
        if (out.flgIssue)                
            flg_skipNow(idxCorecs) = 1;
            thisSkipReason = 'readsac reported issue';
            skipReasonList = [skipReasonList; thisSkipReason];
            fprintf(1,[traceName, '\n\tskipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
            %tripletIdx = 0;
        end
    else
        tripletIdx=0;
    end
    
    
    % -------------------
    % B1. Load waveforms
    % -------------------
    if (~flg_skipNow(itr))
        
        
        % 1. Read current trace
        out = read_sac_trace3(traceFullName);
        if isSocalTrace; s = out.sraw/100;    % Convert to SI units for both acc and vel
        else             s = out.sraw;        % Socal-amps are in [cm/x], Nocal-amps are in [m/x]
        end
        t  = out.t;           
        ns = numel(s);
        sr = out.sr;
        dt = 1/out.sr;
        ts = t(1);
        %[h1] = plot_wform(s,t,[],[],[100 150],591);
        
        TraceList.station.sr(itr) = sr;
        TraceList.eq.ts(itr)      = ts;
        
        if rem(sr,2)~=0; fnameList = [fnameList; traceFullName]; end
        
        % 2. Read vertical component (for picking)
        if tripletIdx==1 
            vertFullName = TraceList.fullName{vertIdx};
            [outpx]      = read_sac_trace3(vertFullName);
            if isSocalTrace; srawpx = outpx.sraw/100;    
            else             srawpx = outpx.sraw;
            end
            tz           = outpx.t;
            nsz          = numel(srawpx);
            srpx         = outpx.sr;
            tzs          = tz(1);
        end
        %clf; hold on; plot(t,s,'w'); plot(tz,srawpx,'c')
        if strcmp(traceName,'72282711.BK.CVS.HNZ.00.sac')
            1+1;
        end

        %fprintf(1,[('sr=',num2str(sr),' Hz, m=',num2str(m),' @hD=',num2str(hD,'%4.0f'),'km)\t'])
        fprintf(1,sprintf(' (sr=%5.1fHz, m=%4.2f @hD=%4.0fkm)\t',sr,m,hD))
        
        % Compute approximate pick location (used by SBPx.m to avoid triggering on foreshocks or pre-event spikes)
        tppx_pred = 0+hD/vp;
        ppxIdx0   = round((tppx_pred-tzs)*srpx)-prePxWindow*srpx;               % Rough guess of ppxIdx minus a couple of seconds
        if ppxIdx0<1; ppxIdx0=1; end
        tAfterPx  = tz(end)-tppx_pred;
        
        % Quality control 2
        % -----------------
        nval   = numel(unique(out.sraw));   % No. of different values. Is too low in some weird waveforms ...
        length = numel(t(t>0))/sr;          % Signal length in seconds after t=0
        
        [~,n_mostFreqVal,~] = mode(out.sraw);   % How often does most frequently occuring amplitude value occur? 
        TraceList.var.v4{itr} = n_mostFreqVal;
        
        if (nval<200)               
            %flg_skipNow(idxCorecs) = 1;
            %flg_skipAfterProc(idxCorecs) = 1;
            thisSkipReason = 'too many identical amps';
            skipReasonList = [skipReasonList; thisSkipReason];
            %TraceList.comment(idxCorecs) = {skipReason};
            fprintf(1,['NOT skipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
        end
        if (length>600)                % ... if longer than 10min
            flg_skipNow(idxCorecs) = 1;
            thisSkipReason              = 'waveform too long';
            skipReasonList = [skipReasonList; thisSkipReason];
            fprintf(1,['skipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
            %tripletIdx = 0;
        end
        if ( (length<20*snpLength) ||(tAfterPx<3) )
            flg_skipNow(idxCorecs) = 1;
            thisSkipReason              = 'waveform too short I';
            skipReasonList = [skipReasonList; thisSkipReason];
            fprintf(1,['skipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
            %tripletIdx = 0;
        end
        if (out.flgIssue ||outpx.flgIssue) 
            flg_skipNow(idxCorecs) = 1;
            thisSkipReason              = 'readsac reported issue';
            skipReasonList = [skipReasonList; thisSkipReason];
            fprintf(1,[traceName, '\n\tskipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
            %tripletIdx = 0;
        end
        if (abs((ns-nsz)/ns)>=0.1) 
            flg_skipNow(idxCorecs) = 1;
            thisSkipReason              = 'Z and H traces have very different lengths';
            skipReasonList = [skipReasonList; thisSkipReason];
            fprintf(1,[traceName, '\n\tskipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
            %tripletIdx = 0;
        end
        if (max(t) < 5)      % ... if time vector is mostly negative
            flg_skipNow(idxCorecs) = 1;
            thisSkipReason              = 'most time-values negative';
            skipReasonList = [skipReasonList; thisSkipReason];
            fprintf(1,['skipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
            %tripletIdx = 0;
        end
        if (sr <= minSr)
            flg_skipNow(idxCorecs) = 1;
            %flg_skipAfterProc(idxCorecs) = 1;
            thisSkipReason                  = 'low sr';
            skipReasonList = [skipReasonList; thisSkipReason];
            fprintf(1,['\n\tskipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
            %tripletIdx = 0;
        end
        if (rem(sr,1) ~= 0)
            flg_skipNow(idxCorecs) = 1;
            thisSkipReason              = 'sr is not integer';
            skipReasonList = [skipReasonList; thisSkipReason];
            fprintf(1,['\n\tskipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
            %tripletIdx = 0;
        end
        if (srpx~=sr)  % Some traces give weird header values with rsac.m, such as sr=6 --> skip'em
            flg_skipNow(idxCorecs) = 1;
            thisSkipReason              = 'more than one sr';
            skipReasonList = [skipReasonList; thisSkipReason];
            fprintf(1,['\n\tskipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
            %tripletIdx = 0;
        end
    else
        tripletIdx=0;
    end
    
     
    % -----------------------------------------------
    % B2. Pick waveform, get absolute pick time, tppx
    % -----------------------------------------------
    % If traces are ok, and absolute pick time is not yet available, get it
    if ( (~flg_skipNow(itr)) && (tripletIdx==1) )
        
        flg_snrIssue = false;
        
        % Prepare vertical trace for picking
        smpx    = srawpx - mean(srawpx(1:2*ntap));                              % Remove mean
        stappx  = taper(smpx,srpx,ntap);                                        % Taper
        if ~isSM; [sh_px]   = bworth(stappx,srpx,fLow_prefilt ,'high',2,'causal');          % Pre-filter (high-pass)
                  [apx,~,~] = bb2accVelDsp(sh_px,srpx,ppxIdx0,'allWform',fLow_prefilt,2,fMode);
                  abpx      = bworth(apx   ,srpx,[fLow_px fUp_px],'band',2,'causal');       
        else      abpx      = bworth(stappx,srpx,[fLow_px fUp_px],'band',2,'causal');       
        end

        % If there is a scsnPx, read it ... 
        if hasScsnPx(itr)
            tppx                     = TraceList.px.p.t(itr);                                 % Absolute time of pick
            [~,ppxIdx_z]             = min(abs(tppx-tz));                                   % What's the closest index to this time on vert. comp.?
            [snr_z,~,~,flg_snrIssue] = get_snr(abpx,ppxIdx_z,srpx,px.Param.signalWindow,px.Param.gapWindow,px.Param.noiseWindow);  % Measure noise before pick on vert. comp.
            if o.plotPx; [h1]        = plot_wform(abpx,tz,tppx,[],[-2 10],591); end
            fprintf(1,[' snr_z: ',num2str(snr_z,'%4.0f'),', '])
            
        % ... if there is not, run sta/lta picker
        else
            
            % ... check if it triggers the STA/LTA pick
            [ppxIdx_z,snr_z,~,~] = SBPx(abpx,1/srpx,ppxIdx0,px.Param,px.Weight,px.Opt);
            fprintf(1,[' snr_z_1: ',num2str(snr_z,'%4.0f'),' '])
            
            % ... 2nd chance, with higher threshold and shift-factor
            if (snr_z<snr_min)
                [ppxIdx_z2,snr_z2,~,~] = SBPx(abpx,1/srpx,ppxIdx0,px.Param,px.Weight,px.Opt);
                fprintf(1,[', snr_z_2: ',num2str(snr_z2,'%4.0f'),' '])
                
                % Keep and use the better one of the two picks
                if snr_z2>snr_z
                    snr_z    = snr_z2;
                    ppxIdx_z = ppxIdx_z2;
                end
            end
        end
        
        % Check if vertical trace has an entry in pxList. If so, assign it to
        % all corecs and overwrite pick fields
        idxPxList = find(cellfun(@(x) ~isempty(x), regexp(vertFullName,pxList.fullName)));
        if ~isempty(idxPxList)
            fprintf(1,'\nOverwriting automatic pick with entry in pxList\n')
            if numel(idxPxList)~=1; fprintf(1,'8ung: more than one pxList-entry found. --> trouble... check.\n'); pause; end
            ppxIdx_z                     = pxList.ppxIdx(idxPxList);
            TraceList.comment(idxCorecs) = {pxList.comment{idxPxList}};
        end

        if (~isempty(ppxIdx_z))
            tppx                            = tz(ppxIdx_z);
            TraceList.px.p.t(idxCorecs)     = tppx;
            TraceList.px.p.hasPx(idxCorecs) = true;
        else
            flg_skipNow(idxCorecs) = 1;
            thisSkipReason         = 'no StaLtaPx'; % (for traces with low SNR, reason will be overwritten)
            skipReasonList         = [skipReasonList; thisSkipReason];
            tppx                   = 0;
            fprintf(1,['\n\tskipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
            %tripletIdx = 0;
        end
        
        % Quality control 3
        % -----------------
        % When the scsn network px makes no sense, snr_z is returned
        % empty from get_noise.m. Set it to zero so that it will be skipped
        if flg_snrIssue; 
            flg_skipNow(idxCorecs)      = 1;
            thisSkipReason              = 'not enough data for SNR';
            skipReasonList              = [skipReasonList; thisSkipReason];
            TraceList.px.p.idx(idxCorecs) = 0;      % Delete pick
            fprintf(1,['\n\tskipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
            %tripletIdx = 0;
        end
        if ( (snr_z>1e11) & ~flg_snrIssue )
            flg_skipNow(idxCorecs)        = 1;
            thisSkipReason                = 'excessive SNR';
            skipReasonList                = [skipReasonList; thisSkipReason];
            TraceList.px.p.idx(idxCorecs) = 0;      % Delete pick
            TraceList.noise.snr(vertIdx)  = snr_z;
            fprintf(1,['\n\tskipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
            %tripletIdx = 0;
        end
    end
    
    
    % -------------------------------------------------------
    % Find pick index on the time vector of actual component 
    % -------------------------------------------------------
    if (~flg_skipNow(itr))
        [~,ppxIdx]              = min(abs(tppx-t));  % Pick index on the time vector of the actual 
        TraceList.px.p.idx(itr) = ppxIdx;            % component (as opposed to ppxIdx_z which is
                                                     % the index the pick has on the time vector of 
                                                     % the verical component
        % Quality control 4
        % -----------------
        % Negative or weird pick times?
        dtpx = tppx-tppx_pred;
        if (ppxIdx<10 ||tppx<t(1) ||tppx<0 ||dtpx>2.5 ||dtpx<-1.5)
            flg_skipNow(idxCorecs)      = 1;
            thisSkipReason              = 'wrong pick';
            skipReasonList              = [skipReasonList; thisSkipReason];
            fprintf(1,['\n\tskipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
            TraceList.px.p.idx(idxCorecs) = 0;          % Delete pick
            %tripletIdx = 0;
        end
        if ((t(end)-tppx)<4*snpLength)               % This one is probably redundant
            flg_skipNow(idxCorecs)      = 1;
            thisSkipReason              = 'waveform too short II';
            skipReasonList              = [skipReasonList; thisSkipReason];
            fprintf(1,['\n\tskipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
            TraceList.px.p.idx(idxCorecs) = 0;      % Delete pick
            %tripletIdx = 0;
        end
        
        
        % Find indices for beginning and end of phases
        % --------------------------------------------
        % Theoretical s-phase arrival
        dtsp                    = hD*sdelay;                    % ts - tp [s]
        tspx                    = tppx + dtsp;
        spxIdx                  = ppxIdx + round(dtsp*sr);
        TraceList.px.s.t(itr)   = single(tspx);
        TraceList.px.s.idx(itr) = spxIdx;
        
        % P-WAVE WINDOW ..............................
        pendIdx = spxIdx-ceil(tpsgap*sr);
        pint    = ppxIdx+1:pendIdx;    % Samples of p-window in full signal ...
        pint    = pint(pint<ns);
        if isempty(pint)
            flg_skipNow(idxCorecs) = 1;
            thisSkipReason         = 'no samples after ppxIdx';
            skipReasonList         = [skipReasonList; thisSkipReason];
            fprintf(1,[traceName, '\n\tskipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
        end

        % S-WAVE WINDOW ..............................
        ns_path = ceil(sr*bt2015_path_duration(hD));           % End of signal according to 
        ns_src  = ceil(sr*get_srcDuration_HanksThatcher72(m)); % conservative (3x) estimates of 
        sendIdx = spxIdx + 3*(ns_path+ns_src);                 % source and path durations
        sint    = pendIdx+1:sendIdx;
        sint    = sint(sint<ns);
        if isempty(sint)
            flg_skipNow(idxCorecs) = 1;
            thisSkipReason         = 'no samples after spxIdx';
            skipReasonList         = [skipReasonList; thisSkipReason];
            fprintf(1,[traceName, '\n\tskipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
        end
        
    else
        tripletIdx=0;
    end
    
    
    % ----------------------------------------------------------
    % BX. Compute acc-, vel- & dsp-waveforms; combine components
    % ----------------------------------------------------------
    % 1. Individual components, highpassed, for fb-amax, PGx, etc.         acc  /vel  /dsp
    % 2. Individual components, bandpassed, for pa, pv & pd                acc2 /vel2 /dsp2
    % 3. Individual components, interpolated to common time vector, for    accti/velti/dspti
    %    envelopes & for saving onset timeseries; truncated 
    % 4. Combined components                                               acc_zen/vel_zen/dsp_zen
    %                                                                      acc_en /vel_en /dsp_en
    if ~flg_skipNow(itr)
        
        % 1. & 2. High- & Bandpassed individual components  ......................................
        % Remove early mean, taper and pre-filter waveform (2-pole high pass)
        preInt = (1:ppxIdx-ntap);
        if isempty(preInt); preInt=(1:ppxIdx); end
        sm   = s - mean(s(preInt));
        %sm   = s - mean(s(1:ppxIdx));
        stap = taper(sm,sr,ntap);
        sh   = bworth(stap,sr,fLow_prefilt     ,'high',fOrder,fMode);
        sb   = bworth(stap,sr,[fLow_prefilt,30],'band',fOrder,fMode); 
        %ns   = numel(sb);
        
        % Integrate to velocity and displacement, and/or differentiate to displacement
        if isSM; [acc ,vel ,dsp ] = sm2accVelDsp(sh,sr,ppxIdx,intMode,fLow_prefilt,fOrder,fMode);
                 [acc2,vel2,dsp2] = sm2accVelDsp(sb,sr,ppxIdx,intMode,fLow_prefilt,fOrder,fMode);
        else     [acc ,vel ,dsp ] = bb2accVelDsp(sh,sr,ppxIdx,intMode,fLow_prefilt,fOrder,fMode);
                 [acc2,vel2,dsp2] = bb2accVelDsp(sb,sr,ppxIdx,intMode,fLow_prefilt,fOrder,fMode);
        end
        mmitmp = pgx2immi(abs(acc),abs(vel));  % Not tested
        mmi    = mmitmp.immi;
        %clf; plot(acc); hold on; plot(ppxIdx,0,'xr'); plot(sendIdx,0,'xr')
        
        
        % 3. Interpolated individual components  .................................................
        % Get envelopes and acc-onset with 100sps, starting from 1 sample
        % before pick; "t" stands for truncated; "i" stands for interpolated
        %interval = (t>=tppx-1/sr &t<(tppx+tmax-1/sr));
        allint   = pint(1):sint(end); % beginning of p-phase to end of s-phase
        tt       = t  (allint);
        acct     = acc(allint);
        velt     = vel(allint);
        dspt     = dsp(allint);
        rawt     = sm (allint);
        % clf; hold on; plot(t,acc2,'k'); plot(tt,acct,'r','lineWidth',2); plot(tt,acct,':k')
        
        % Interpolate to 100sps
        tti   = tt(1):.01:tt(end);
        accti = interp1(tt,acct,tti);
        velti = interp1(tt,velt,tti);
        dspti = interp1(tt,dspt,tti);
        % hold on; plot(tti,accti,'dm')

        
        % 4. Combined components  ...............................................................
        if tripletIdx==1; accMat = cell(3,1);
                          velMat = cell(3,1);
                          dspMat = cell(3,1);
        end
        accMat{tripletIdx} = abs(acc(allint));
        velMat{tripletIdx} = abs(vel(allint));
        dspMat{tripletIdx} = abs(dsp(allint));
        %clf; plot(pint,0,'xr'); hold on; plot(sint,0,'xb'); plot(ppxIdx,0,'db'); plot(spxIdx,0,'dr'); plot(acc,'k')
        
        if tripletIdx==3;
            
            % Which one is vertical component?
            iiH = setdiff(1:3,zHit);   % triplet indices of the two horizonal traces 
            if ~strcmp(TraceList.station.ocode(idxCorecs(zHit)),'Z');
                fprintf(1,'8UNG: does zHit really tell you which of the triplet is the vertical trace?\n ')
                cvert = cvert+1;
            end
            
            % Z/E/N waveforms may not have equal lengths; find shortest one and truncate others to that length
            nsminAcc = min(cellfun(@(x) size(x,1), accMat));
            nsminVel = min(cellfun(@(x) size(x,1), velMat));
            nsminDsp = min(cellfun(@(x) size(x,1), dspMat));
            nsmin3   = min([nsminAcc,nsminVel,nsminDsp]);
            accMat2  = cell2mat(cellfun(@(x) x(1:nsmin3),accMat','uniformOutput',0));
            velMat2  = cell2mat(cellfun(@(x) x(1:nsmin3),velMat','uniformOutput',0));
            dspMat2  = cell2mat(cellfun(@(x) x(1:nsmin3),dspMat','uniformOutput',0));
            
            acc_zen = sqrt(sum(accMat2.^2,2));
            vel_zen = sqrt(sum(velMat2.^2,2));
            dsp_zen = sqrt(sum(dspMat2.^2,2));
            mmitmp  = pgx2immi(abs(acc_zen),abs(vel_zen));
            mmi_zen = mmitmp.immi;
            
            acc_en  = sqrt(sum(accMat2(:,iiH).^2,2));
            vel_en  = sqrt(sum(velMat2(:,iiH).^2,2));
            dsp_en  = sqrt(sum(dspMat2(:,iiH).^2,2));
            mmitmp  = pgx2immi(abs(acc_en),abs(vel_en));
            mmi_en  = mmitmp.immi;
            %clf; hold on; plot(t,abs(acc),'k'); plot(tt,acc_zen,'r'); plot(tt,acc_en ,'b')
            
            % Windows relative to ppxIdx, rather than relative to signal start
            pintrel   = pint-pint(1)+1;
            sintrel   = sint-pint(1)+1;
            allintrel = pintrel(1):sintrel(end);

            pintrel   = pintrel  (pintrel  <=nsmin3);
            sintrel   = sintrel  (sintrel  <=nsmin3);
            allintrel = allintrel(allintrel<=nsmin3);
        end
    else
        tripletIdx=0;
    end
    
    
    
    % ------------------------------------------------------
    % BX. Process Waveforms / Compute Time-Series Statistics
    % ------------------------------------------------------
    % 1. Final Peak Amplitudes
    % 2. Peak Amplitude Evolution
    % 3. Filter Bank Amplitudes
    % 4. Times to Reach Threshold Amplitudes
    % 5. Envelopes                 
    % 6. Save initial waveforms  
	% 7. Miscellaneous           
    % 8. GM Outlier Check        
    % 9. Write skipReasonList to commment field
    % 10.Various Plots
    if (~flg_skipNow(itr))

        
        % 1. Final Peak Amplitudes ...........................................................
        
        % Measure max amps in P- & S-phase windows separately on individual components
        [ppga,ppgaIdx]            = max(abs(acc(pint)));
        [ppgv,ppgvIdx]            = max(abs(vel(pint)));
        [ppgd,ppgdIdx]            = max(abs(dsp(pint)));
        TraceList.pga.p.ampi(itr) = ppga;
        TraceList.pgv.p.ampi(itr) = ppgv;
        TraceList.pgd.p.ampi(itr) = ppgd;
        TraceList.pga.p.idxi(itr) = ppgaIdx;
        TraceList.pgv.p.idxi(itr) = ppgvIdx;
        TraceList.pgd.p.idxi(itr) = ppgdIdx;
        
        
        % Make S-window indices rel. to ppxIdx. Add 1 because its a fencepost problem, 
        dsIdx = (pendIdx-ppxIdx);  % and 1 because the S-window starts at pendIdx+1.
        
        [spga,spgaIdx]            = max(abs(acc(sint)));
        [spgv,spgvIdx]            = max(abs(vel(sint)));
        [spgd,spgdIdx]            = max(abs(dsp(sint)));
        spgaIdx                   = spgaIdx+dsIdx;
        spgvIdx                   = spgvIdx+dsIdx;
        spgdIdx                   = spgdIdx+dsIdx;
        TraceList.pga.s.ampi(itr) = spga;
        TraceList.pgv.s.ampi(itr) = spgv;
        TraceList.pgd.s.ampi(itr) = spgd;
        TraceList.pga.s.idxi(itr) = spgaIdx;
        TraceList.pgv.s.idxi(itr) = spgvIdx;
        TraceList.pgd.s.idxi(itr) = spgdIdx;
        %clf; hold on;  plot(acc,'w'); plot(ppxIdx+pgaIdx,acc(ppxIdx+pgaIdx),'or');
        
        if spga>=ppga; TraceList.pga.pns.ampi(itr) = spga;
                       TraceList.pga.pns.idxi(itr) = spgaIdx;
        else           TraceList.pga.pns.ampi(itr) = ppga;
                       TraceList.pga.pns.idxi(itr) = ppgaIdx;
        end
        if spgv>=ppgv; TraceList.pgv.pns.ampi(itr) = spgv;
                       TraceList.pgv.pns.idxi(itr) = spgvIdx;
        else           TraceList.pgv.pns.ampi(itr) = ppgv;
                       TraceList.pgv.pns.idxi(itr) = ppgvIdx;
        end
        if spgd>=ppgd; TraceList.pgd.pns.ampi(itr) = spgd;
                       TraceList.pgd.pns.idxi(itr) = spgdIdx;
        else           TraceList.pgd.pns.ampi(itr) = ppgd;
                       TraceList.pgd.pns.idxi(itr) = ppgdIdx;
        end
        pgv = TraceList.pgv.pns.ampi(itr); % Used for outlier check
        
         
        % Measure max amps in S-wave window  on combined components...
        if tripletIdx==3;
            
            % ... all three components ....................................
            [ppga,ppgaIdx]                  = max(abs(acc_zen(pintrel)));
            [ppgv,ppgvIdx]                  = max(abs(vel_zen(pintrel)));
            [ppgd,ppgdIdx]                  = max(abs(dsp_zen(pintrel)));
            TraceList.pga.p.ampzen(vertIdx) = ppga;
            TraceList.pgv.p.ampzen(vertIdx) = ppgv;
            TraceList.pgd.p.ampzen(vertIdx) = ppgd;
            TraceList.pga.p.idxzen(vertIdx) = ppgaIdx;
            TraceList.pgv.p.idxzen(vertIdx) = ppgvIdx;
            TraceList.pgd.p.idxzen(vertIdx) = ppgdIdx;
            
            % Make S-window indices rel. to ppxIdx. Add 1 because its a fencepost problem,
            dsIdx = (pendIdx-ppxIdx);  % and 1 because the S-window starts at pendIdx+1.

            [spga,spgaIdx]                  = max(abs(acc_zen(sintrel)));
            [spgv,spgvIdx]                  = max(abs(vel_zen(sintrel)));
            [spgd,spgdIdx]                  = max(abs(dsp_zen(sintrel)));
            spgaIdx                         = spgaIdx+dsIdx;
            spgvIdx                         = spgvIdx+dsIdx;
            spgdIdx                         = spgdIdx+dsIdx;
            TraceList.pga.s.ampzen(vertIdx) = spga;
            TraceList.pgv.s.ampzen(vertIdx) = spgv;
            TraceList.pgd.s.ampzen(vertIdx) = spgd;
            TraceList.pga.s.idxzen(vertIdx) = spgaIdx;
            TraceList.pgv.s.idxzen(vertIdx) = spgvIdx;
            TraceList.pgd.s.idxzen(vertIdx) = spgdIdx;
            
            if spga>=ppga; TraceList.pga.pns.ampzen(vertIdx) = spga;
                           TraceList.pga.pns.idxzen(vertIdx) = spgaIdx;
            else           TraceList.pga.pns.ampzen(vertIdx) = ppga;
                           TraceList.pga.pns.idxzen(vertIdx) = ppgaIdx;
            end
            if spgv>=ppgv; TraceList.pgv.pns.ampzen(vertIdx) = spgv;
                           TraceList.pgv.pns.idxzen(vertIdx) = spgvIdx;
            else           TraceList.pgv.pns.ampzen(vertIdx) = ppgv;
                           TraceList.pgv.pns.idxzen(vertIdx) = ppgvIdx;
            end
            if spgd>=ppgd; TraceList.pgd.pns.ampzen(vertIdx) = spgd;
                           TraceList.pgd.pns.idxzen(vertIdx) = spgdIdx;
            else           TraceList.pgd.pns.ampzen(vertIdx) = ppgd;
                           TraceList.pgd.pns.idxzen(vertIdx) = ppgdIdx;
            end
            TraceList.mmi.pns.ampzen(vertIdx) = max(mmi_zen);
            
            
            % ... horizontal components ...................................
            [ppga,ppgaIdx]                 = max(abs(acc_en(pintrel)));
            [ppgv,ppgvIdx]                 = max(abs(vel_en(pintrel)));
            [ppgd,ppgdIdx]                 = max(abs(dsp_en(pintrel)));
            TraceList.pga.p.ampen(vertIdx) = ppga;
            TraceList.pgv.p.ampen(vertIdx) = ppgv;
            TraceList.pgd.p.ampen(vertIdx) = ppgd;
            TraceList.pga.p.idxen(vertIdx) = ppgaIdx;
            TraceList.pgv.p.idxen(vertIdx) = ppgvIdx;
            TraceList.pgd.p.idxen(vertIdx) = ppgdIdx;
            
            [spga,spgaIdx]                 = max(abs(acc_en(sintrel)));
            [spgv,spgvIdx]                 = max(abs(vel_en(sintrel)));
            [spgd,spgdIdx]                 = max(abs(dsp_en(sintrel)));
            spgaIdx                        = spgaIdx+dsIdx;
            spgvIdx                        = spgvIdx+dsIdx;
            spgdIdx                        = spgdIdx+dsIdx;
            TraceList.pga.s.ampen(vertIdx) = spga;
            TraceList.pgv.s.ampen(vertIdx) = spgv;
            TraceList.pgd.s.ampen(vertIdx) = spgd;
            TraceList.pga.s.idxen(vertIdx) = spgaIdx;
            TraceList.pgv.s.idxen(vertIdx) = spgvIdx;
            TraceList.pgd.s.idxen(vertIdx) = spgdIdx;
            
            if spga>=ppga; TraceList.pga.pns.ampen(vertIdx) = spga;
                           TraceList.pga.pns.idxen(vertIdx) = spgaIdx;
            else           TraceList.pga.pns.ampen(vertIdx) = ppga;
                           TraceList.pga.pns.idxen(vertIdx) = ppgaIdx;
            end
            if spgv>=ppgv; TraceList.pgv.pns.ampen(vertIdx) = spgv;
                           TraceList.pgv.pns.idxen(vertIdx) = spgvIdx;
            else           TraceList.pgv.pns.ampen(vertIdx) = ppgv;
                           TraceList.pgv.pns.idxen(vertIdx) = ppgvIdx;
            end
            if spgd>=ppgd; TraceList.pgd.pns.ampen(vertIdx) = spgd;
                           TraceList.pgd.pns.idxen(vertIdx) = spgdIdx;
            else           TraceList.pgd.pns.ampen(vertIdx) = ppgd;
                           TraceList.pgd.pns.idxen(vertIdx) = ppgdIdx;
            end
            TraceList.mmi.pns.ampen(vertIdx) = max(mmi_en);
        end
        
        
        % 2. Peak Amplitude Evolution ...........................................................
        if o.useBP; [pa,pna,paIdx,~,~] = measure_tdpa(acc2,ppxIdx,sr,tnoise,tmax,snpLength,[]);
                    [pv,pnv,pvIdx,~,~] = measure_tdpa(vel2,ppxIdx,sr,tnoise,tmax,snpLength,[]);
                    [pd,pnd,pdIdx,~,~] = measure_tdpa(dsp2,ppxIdx,sr,tnoise,tmax,snpLength,[]);
        else        [pa,pna,paIdx,~,~] = measure_tdpa(acc ,ppxIdx,sr,tnoise,tmax,snpLength,[]);
                    [pv,pnv,pvIdx,~,~] = measure_tdpa(vel ,ppxIdx,sr,tnoise,tmax,snpLength,[]);
                    [pd,pnd,pdIdx,~,~] = measure_tdpa(dsp ,ppxIdx,sr,tnoise,tmax,snpLength,[]);
        end
        %plot_wform_and_tdpa(acc2,t,tppx,[],[tppx-1.1*tnoise 1.1*tppx+tmax],pa,snpLength,1115)
        pitmp                  = pgx2immi(pa,pv);
        pi                     = trim_amax(pitmp.immi,nsnpmin);
        pa                     = trim_amax(pa,nsnpmin);
        pv                     = trim_amax(pv,nsnpmin);
        pd                     = trim_amax(pd,nsnpmin);
        TraceList.pga.tsi{itr} = single(pa);
        TraceList.pgv.tsi{itr} = single(pv);
        TraceList.pgd.tsi{itr} = single(pd);
        TraceList.mmi.tsi{itr} = single(pi);
        
        % Measure pa/pv/pd on combined components...
        if tripletIdx==3;
            [pa,pna,paIdx,~,~]           = measure_tdpa(acc_zen,1,sr,tnoise,tmax,snpLength,[]);
            [pv,pnv,pvIdx,~,~]           = measure_tdpa(vel_zen,1,sr,tnoise,tmax,snpLength,[]);
            [pd,pnd,pdIdx,~,~]           = measure_tdpa(dsp_zen,1,sr,tnoise,tmax,snpLength,[]);
            pitmp                        = pgx2immi(pa,pv);
            pi                           = trim_amax(pitmp.immi,nsnpmin);
            pa                           = trim_amax(pa,nsnpmin);
            pv                           = trim_amax(pv,nsnpmin);
            pd                           = trim_amax(pd,nsnpmin);
            TraceList.pga.tszen{vertIdx} = single(pa);
            TraceList.pgv.tszen{vertIdx} = single(pv);
            TraceList.pgd.tszen{vertIdx} = single(pd);
            TraceList.mmi.tszen{vertIdx} = single(pi);
        end
        
        
        % 3. Filter Bank Amplitudes .............................................................
        % pass velocity waveform through | BUTTERWORTH FILTER BANK |
        if ~o.omitAmax
            [amax,~,amaxIdx,sout,~] = measure_tdpa(vel,ppxIdx,sr,0,tmax,snpLength,fc);
            amax                    = trim_amax(amax,nsnpmin);
            TraceList.fb.amax{itr}  = single(amax);
            %TraceList.nbnoise{itr} = single(nbNoise);
        end
        
        
        % X. Azimuth and waveform polarisation  .................................................
        if o.estimateBaz
            if tripletIdx==3
                
                % Different methods are somewhat described in m_tools/azimuth/readme
                
                % True backazimuth  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                %az  = azimuth(eqLat,eqLon,stLat,stLon);
                bAz          = azimuth(stLat,stLon,eqLat,eqLon);
                bazi.trueVal = bAz;
                %[bAz,ray] = get_backAzimuth(eqLat,eqLon,stLat,stLon);
                
                
                % M1: Method implemented in onsite framework  . . . . . . . . . . . . . . . . . .
                Z = dspMat2(:,zHit)'; % TraceList.station.ocode(idxCorecs)
                E = dspMat2(:,eHit)';
                N = dspMat2(:,nHit)';
                
                bazi.onsite.baziHat = zeros(nbaz,1);
                for iint = 1:nbaz
                    int                       = int8(1:bazIntervalVector(iint)*sr);
                    bazi.onsite.baziHat(iint) = dsp2bazi_onsite(Z,E,N,int);
                end
                
                % M1b: Polarization with Zach Ross' method
                
                
                % M2: Sac method  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                bazi.sac.baziHat = zeros(nbaz,1);
                bazi.sac.incHat = zeros(nbaz,1);
                bazi.sac.polHat = zeros(nbaz,1);
                for iint = 1:nbaz
                    int                 = int8(1:bazIntervalVector(iint)*sr);
                    [az,incid,polariz]  = getWinPolariz_MAM(Z,E,N,int);
                    
                    backazimuth = az-180;
                    if backazimuth<0; backazimuth=backazimuth+360; end
                    
                    bazi.sac.baziHat(iint) = backazimuth;
                    bazi.sac.incHat(iint)  = incid;
                    bazi.sac.polHat(iint)  = polariz;
                end
                
                
                % M3: Nicolas Deichmanns method   . . . . . . . . . . . . . . . . . . . . . . . .
                bazi.deichmann.baziHat = zeros(nbaz,1);
                bazi.deichmann.incHat  = zeros(nbaz,1);
                for iint = 1:nbaz
                    int      = int8(1:bazIntervalVector(iint)*sr);
                    [az,inc] = azimuth_incidence_ND(Z,E,N,int);
                    
                    backazimuth = az-180;
                    if backazimuth<0; backazimuth=backazimuth+360; end
                    
                    bazi.deichmann.baziHat(iint) = backazimuth;
                    bazi.deichmann.incHat(iint)  = inc;
                end
                
                
                % M4: Lockman and Allen, 2005, BSSA   . . . . . . . . . . . . . . . . . . . . . .
                bazi.la05.baziHat = zeros(nbaz,1);
                for iint = 1:nbaz
                    int                    = int8(1:bazIntervalVector(iint)*sr);
                    az                     = dsp2azi_la05(Z,E,N,int,0.5);
                    backazimuth = az-180;
                    if backazimuth<0; backazimuth=backazimuth+360; end
                    bazi.la05.baziHat(iint) = az;
                end
                
                
                % M5: Method used by seismicHandler   . . . . . . . . . . . . . . . . . . . . . .
                % Other methods ...   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                TraceList.scalFeature{vertIdx}.bazi = bazi; 
            end
        end
        
        % X. P/S discriminant  ..................................................................
        
        
        % 4. Times to Reach Threshold Amplitudes  ...............................................
        % ... on individual components
        idxAcc     = get_threshold_indices(acc(ppxIdx:end),accThresholds );
        idxVel     = get_threshold_indices(vel(ppxIdx:end),velThresholds );
        idxImmi    = get_threshold_indices(mmi(ppxIdx:end),immiThresholds);
        tOcc.acc1  = idxAcc;
        tOcc.vel1  = idxVel;
        tOcc.mmi1  = idxImmi;
        
        TraceList.scalFeature{itr}.tOcc = tOcc;

        if tripletIdx==3;
        
            tOcc = TraceList.scalFeature{vertIdx}.tOcc;

            % ... on vector sum of all three components
            idxAcc       = get_threshold_indices(acc_zen(1:end),accThresholds);
            idxVel       = get_threshold_indices(vel_zen(1:end),velThresholds);
            idxImmi      = get_threshold_indices(mmi_zen(1:end),immiThresholds);
            tOcc.acc_zen = idxAcc;
            tOcc.vel_zen = idxVel;
            tOcc.mmi_zen = idxImmi;

            % ... on vector sum of horizontal components
            idxAcc      = get_threshold_indices(acc_en(1:end),accThresholds);
            idxVel      = get_threshold_indices(vel_en(1:end),velThresholds);
            idxImmi     = get_threshold_indices(mmi_en(1:end),immiThresholds);
            tOcc.acc_en = idxAcc;
            tOcc.vel_en = idxVel;
            tOcc.mmi_en = idxImmi;
            
            TraceList.scalFeature{vertIdx}.tOcc = tOcc;
        end
        
        
        % X. Scalar waveform features .....................................
        if o.compute_scalFeatures
            
            % Use dsp-field for raw time series
            [iFeatures,varout] = get_vsumFeatures_growingSnippet(acct,velt,rawt,sr,snpLength,nsnpFeatures);
            TraceList.scalFeature{itr}.acci.zcr  = iFeatures.acc.zcr;
            TraceList.scalFeature{itr}.acci.skew = iFeatures.acc.skewness;
            TraceList.scalFeature{itr}.acci.kurt = iFeatures.acc.kurtosis;
            TraceList.scalFeature{itr}.acci.qtr  = iFeatures.acc.qtr;
            TraceList.scalFeature{itr}.acci.madV = iFeatures.acc.mad;
            TraceList.scalFeature{itr}.acci.f38V = iFeatures.acc.f38;
            TraceList.scalFeature{itr}.acci.k2V  = iFeatures.acc.k2;
            
            TraceList.scalFeature{itr}.veli.cav  = iFeatures.vel.cax;
            TraceList.scalFeature{itr}.rawi.zcr  = iFeatures.dsp.zcr;
            TraceList.scalFeature{itr}.rawi.skew = iFeatures.dsp.skewness;
            TraceList.scalFeature{itr}.rawi.kurt = iFeatures.dsp.kurtosis;
            TraceList.scalFeature{itr}.rawi.qtr  = iFeatures.dsp.qtr;
            TraceList.scalFeature{itr}.rawi.cav  = iFeatures.dsp.cax;
            
            taucOut = get_tauC(velt,dspt);
            tauC    = taucOut.tau_c(varout.iE);
            TraceList.scalFeature{itr}.veli.tauc = tauC;
            
            int1 = 1:0.2/dt;
            int2 = (0.2+dt)/dt:0.4/dt;
            rvar = var(acc(int2))/var(acc(int1));
            TraceList.scalFeature{itr}.acci.rvar = rvar;
            
            if tripletIdx==3
                
                % Feature Z/H: max(z)/max(en)
                accz = accMat{zHit}; %clf; plot(accz,'k');hold on;plot(acc_en,'r')
                zhr  = get_zhr_feature(accz,acc_en,sr,snpLength,nsnpFeatures);
                TraceList.scalFeature{vertIdx}.acczen.zhr = zhr;
                
                % Feature maxStep: max([diff(accMat{1}); diff(acc_e); diff(acc_n)]); or on raw?
                maxstep = get_maxstep_feature(rawt,sr,snpLength,nsnpFeatures);
                TraceList.scalFeature{vertIdx}.rawi.maxStep = maxstep.val;
            end
        end

        
        
        % 5. Envelopes  ..........................................................................
        envNsWidth          = round(envWidth*100);
        envNsShift          = round(envShift*100);
        [accEnv,~]          = get_numeric_wform_envelope(accti,envNsWidth,envNsShift);
        [velEnv,~]          = get_numeric_wform_envelope(velti,envNsWidth,envNsShift);
        [dspEnv,~]          = get_numeric_wform_envelope(dspti,envNsWidth,envNsShift);
        TraceList.var.v6{itr} = [single(accEnv)'; single(velEnv)'; single(dspEnv)'];
        %TraceList.accNoise{itr} = single(accEnv);
        % gcf; clf; hold on; plot(t-tppx,abs(acc2),'k') plot(t(idx),accEnv,'r')

        
        % 6. Save initial waveforms  .............................................................
        % Save 100sps acc for first couple of seconds
        interval            = (tti<(tppx+te_onset-1/sr));
        accOnset            = accti(interval);
        velOnset            = velti(interval);
        dspOnset            = dspti(interval);
        TraceList.var.v7{itr} = [single(accOnset);single(velOnset);single(dspOnset)];
        %tOnset   = tti(interval); clf; hold on; plot(tti,accti,'k'); plot(tOnset,accOnset,'r')

        
        % 7. Miscellaneous  .............................................................
        
        % Measure width of 1/4 of first displacement pulse by finding first
        % zero-crossing in velocity trace
        [pulseWidthVel,~]   = get_firstPulseWidth(vel2,ppxIdx,sr);
        [pulseWidthAcc,~]   = get_firstPulseWidth(acc2,ppxIdx,sr);
        TraceList.var.v2{itr} = [pulseWidthVel,pulseWidthAcc];
         
        % Measure and save SNR on bandpassed acceleration trace
        [snr,~,~,~]              = get_snr(acc2,ppxIdx,sr,px.Param.signalWindow,px.Param.gapWindow,px.Param.noiseWindow);
        TraceList.noise.snr(itr) = single(snr);
        
        % Measure and save acceleraton noise amps
        noiseAmps               = acc2(ntap+1:ppxIdx);
        naPctl                  = prctile(noiseAmps,[5, 50, 95]);
        TraceList.noise.acc{itr} = single(naPctl);
        
        if (naPctl(3)==naPctl(2) |abs(naPctl(3))<eps )
            flg_skipNow(idxCorecs) = 1;
            skipReason             = 'something wrong with noise amplitudes';
            skipReasonList = [skipReasonList; thisSkipReason];
            fprintf(1,['skipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
            %tripletIdx = 0;
        end
        % TraceList.dspNoise{itr} = single(pnd);
        % TraceList.velNoise{itr} = single(pnv);
        % TraceList.accNoise{itr} = single(naPctl);
        

        % 8. GM Outlier Check  ..........................................................
        % For large events, pgv may not come from nucleation point
        rld    = 10.^(a+b*m); % Rupture length from Wells and Coppersmith 1994
        hD_min = hD - rld;  if (hD_min<0); hD_min=0; end
        hD_max = hD + rld;
        
        % Compute CH2007 GMPEs for comparison
        if (strcmp(orientationCode,'Z'))
            [~,up,~,~]   = CH2007_nSigma(m,hD_min,'PGV','Z','S','R',0,nsigma);
            [~,~,lo,~]   = CH2007_nSigma(m,hD_max,'PGV','Z','S','R',0,nsigma);
        elseif ( strcmp(orientationCode,'N')||(strcmp(orientationCode,'E')) )
            [~,up,~,~]   = CH2007_nSigma(m,hD_min,'PGV','H','S','R',0,nsigma);
            [~,~,lo,~]   = CH2007_nSigma(m,hD_max,'PGV','H','S','R',0,nsigma);
        else
            fprintf(1,'Something wrong with orientationCode, check. now.\n')
        end
        
        up                    = up/1e2;    % Convert from [cm/s] to [m/s]
        lo                    = lo/1e2;
        TraceList.pgv.up(itr) = up;
        TraceList.pgv.lo(itr) = lo;
        
        % ... if it is, set skip-flag
        if ( (pgv>up) || (pgv<lo) )
            %flg_skipNow(idxCorecs) = 1;
            %flg_skipAfterProc(idxCorecs) = 1;
            thisSkipReason                  = 'outlier';
            skipReasonList = [skipReasonList; thisSkipReason];
            %TraceList.comment(idxCorecs) = {skipReason};
            fprintf(1,['\n\tNOT skipping ',num2str(numel(idxCorecs)),' traces: ',thisSkipReason])
        end
        
        
        % 9. Write skipReasonList to commment field of all corecorded traces  ................
        TraceList.comment(idxCorecs) = {unique(skipReasonList)};

        
        % 10. Various Plots ............................................................
        if o.plotWform &&tripletIdx==3
            idxH = setdiff(idxCorecs,vertIdx);
            out  = plot_wforms_3comps_from_same_trList(vertIdx,idxH(1),idxH(2),TraceList);
        % clf; hold on; plot(out.t3,out.acc_zen,'k','lineWidth',2); plot(t(allint),acc_zen,'r')
        end
        %             % Find indices of plotting window
        %             tppx        = t(ppxIdx);
        %             sIdx        = round(ppxIdx - T*sr*tFrac);      % t-index of window-start
        %             eIdx        = round(ppxIdx + T*sr*(1-tFrac));  % t-index of window-end
        %             titleString = strcat(['Event ',eqName, ' (M',num2str(m),') on ',nw,'.',sta,'.',chan, ...
        %                 ' at distance ',num2str(hD),'km. t0: ',TraceList.eq.date{itr},', ',TraceList.eq.t0{itr}]);
        %             [h1] = plot_nbWforms(vel,sout,amaxIdx,t,tppx,tspx,velNoise1,nbNoise,sIdx,eIdx,fc,titleString,tFrac,T);
        %             [h2] = plot_nbPGV(amax,fc,nbNoise,[lo,pgv,up],snpLength,titleString,o.verbose);
        %             plot_nbPGV_tseries(TraceList.fb.amax{itr},TraceList.nbnoise{itr},fc,snpLength,titleString)
        %             plot_wform(vel,t,tppx,tspx,[-1 20], 592)
        %             %plotTraceSummary(traceFullName,TraceList,8,1,49)
        
        
        
    else
        tripletIdx=0;
        
        % Add reason. These entries will later be relocated to SkipList
        TraceList.comment(idxCorecs) = {unique(skipReasonList)};
        %TraceList.comment{itr} = skipReasonList;
        
        % Measure and save acceleraton noise amps anyways...
        % Notice that these noise amps are not directly comparable to those
        % of the processed (non-skipped) traces since they have differfent
        % filtering/processing
        if exist('s','var')
            iS = ntap;
            iE = ntap+2*round(sr);
            if ~isSM;       s  = diff(s)*sr; end
            if iE>numel(s); iE = numel(s);   end
            noiseAmps = s(iS:iE);   % <s> is the raw, gain-corrected, signal
            naPctl    = prctile(noiseAmps-mean(noiseAmps),[5, 50, 95]);
            TraceList.noise.acc{itr} = single(naPctl);
        end
        
        if (o.plotQDP)
            if ~exist('ppxIdx','var'); ppxIdx=[]; end
            figure(513); clf; hold on
            if (~isempty(ppxIdx)); plot_wform(scounts,tcounts,[],[],592)
            else                   plot_wform(scounts,tcounts,tcounts(ppxIdx),[],592)
            end 
            ylabel('Counts','fontSize',ftSize)
            xlabel('No. of samples','fontSize',ftSize)
            title(['SKIPPED: ',skipReasonList{1}],'fontSize',ftSize)
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            line([xl(1) xl(2)],[yl(1) yl(2)],'color','r','lineStyle',':')
            line([xl(1) xl(2)],[yl(2) yl(1)],'color','r','lineStyle',':')
        end
    end
    
    % If all three traces of triplet have been processed, reset tripletIdx
    if tripletIdx==3; tripletIdx = 0; end
end

% Add all traces on blackList to SkipList
% fprintf(1,'\n\nidentifying blacklisted traces...\n')
% fprintf(1,'\n\nnew blacklist function not yet checked. do it. now.\n')
% pause
% 
% [blackIdxList,~,blComments]       = find_blackListed_traces(TraceList,blackList);
% flg_skipBlackListed               = false(ntr_bak,1);
% flg_skipBlackListed(blackIdxList) = true;
% nbl = numel(blackIdxList);
% for ibl = 1:nbl
%     skipReasonList                       = [TraceList.comment{ibl}; blComments{ibl}];
%     TraceList.comment{blackIdxList(ibl)} = skipReasonList;
% end

% nbl = numel(blackList.eq.m);
% for ibl = 1:nbl
%     print_iteration_numbers(ibl,nbl,'hundreds')
%     blFullName   = blackList.fullName{ibl};
%     idxBlackListed = find(cellfun(@(x) ~isempty(x), regexp(blFullName,TraceList.fullName)));
%     if ~isempty(idxBlackListed)
%     
%         % Identify corecords of black-listed traces
%         traceFullName = TraceList.fullName{idxBlackListed};
%         recordName    = get_recordName(traceFullName);    
%         hits          = regexp(TraceList.fullName,recordName);
%         idxCorecs     = find(cellfun(@(x) ~isempty(x), hits));
%         nCorec        = numel(idxCorecs);
%         if nCorec~=3; fprintf(1,'8ung: not three corecs found. check.\n'); pause; end
%         flg_skipBlackListed(idxCorecs) = true;
%         %TraceList.comment(idxCorecs) = {blackList.comment{ibl}};
%         for ii = 1:numel(idxCorecs)
%             skipReasonList                   = [TraceList.comment{idxCorecs(ii)}; blackList.comment{ibl}];
%             TraceList.comment{idxCorecs(ii)} = skipReasonList;
%         end
%     end
% end
%fprintf(1,sprintf(' done. %i records flagged for skipping because they were black-listed.\n\n',sum(flg_skipBlackListed)))


% Cut out the last 1-2 records if not all it's corecs have been processed,
% e.g. because the no. of traces was limited by ntrMaxAbs
if tripletIdx~=0
    fprintf(1,'Check if last 1-2 traces are ok, or if they lack corec.\n')
    pause
end

flg_skip_tooMany            = false(ntr_bak,1);
flg_skip_tooMany(ntr+1:end) = true;
idx_tooMany                 = find(flg_skip_tooMany);
for id = 1:numel(idx_tooMany)
    skipReasonList                     = [TraceList.comment{idx_tooMany(id)}; 'more traces than ntrMaxAbs'];
    TraceList.comment{idx_tooMany(id)} = skipReasonList;
end
%TraceList.comment(flg_skip_tooMany) = {'more traces than ntrMaxAbs'};

% Add all traces with flg_skipNow(itr) = true to SkipList
flg_skip = (flg_skipNow |flg_skipAfterProc |flg_skip_tooMany);
skipIdx  = find( flg_skip);
keepIdx  = find(~flg_skip);
fprintf(1,['==> adding a total of ',num2str(numel(skipIdx)),' trace to SkipList ...'])
tmpSkipList = TraceList.selectSubList(skipIdx);
SkipList.appendList(tmpSkipList); clear tmpSkipList
TraceList = TraceList.selectSubList(keepIdx);       % Remove skipped traces and outliers from traceList
fprintf(1,['    done.\n'])

% Add a comment for all SM records with an unclipped colocated BB record
fprintf(1,'Finding double records (not skipping) ...')
idx_sm2    = find_double_records(TraceList);
% lgc          = false(numel(TraceList.eq.m),1);
% lgc(skipIdx) = true;
% idx_bb       = find(lgc);
for id = 1:numel(idx_sm2)
    skipReasonList                 = [TraceList.comment{idx_sm2(id)}; 'has unclipped BB record'];
    TraceList.comment{idx_sm2(id)} = skipReasonList;
end

% Make sure that number of traces is multiple of three. THIS SHOULD BE OBSOLETE
ntr_fin = numel(TraceList.eq.m);
ntr_rem = rem(ntr_fin,3);
if ntr_rem~=0
    fprintf(1,'ntr is not multiple of three... why? check. now. do it. And press enter to remove last entry or two...\n')
    pause
    TraceList = TraceList.selectSubList(1:ntr_fin-ntr_rem);
end

%TraceList.comment(lgc) = {'has unclipped BB record'};
%nskp         = numel(skipIdx);
%fprintf(1,['    adding another ',num2str(nskp),' SM traces which have unclipped BB corecs to SkipList ...'])
%tmpSkipList  = TraceList.selectSubList(skipIdx);
%tmpSkipList.comment([1:nskp]') = {'has unclipped BB record'};
%SkipList.appendList(tmpSkipList); clear tmpSkipList
%TraceList.removeSkipped(lgc);
fprintf(1,['    done.\n'])

SkipList.consolidate;

% Print and save summary
summary_trList(TraceList,SkipList);

% No. of cases where tref from sacfile was not equal origin time t0
%c_absTimeIssue = sum(find(cellfun(@(x) ~isempty(x), TraceList.comment)));
fprintf(1,['\n',num2str(ctref),' traces have tref not equal t0\n'])
fprintf(1,[num2str(cvert),' times was there a problem identifying the vertical of the triplet\n'])



% ---------------------------
% B6. Save filter bank output
% ---------------------------
if o.saveOut
    
    fprintf(1,'\n   ###  Storing filter-bank output to harddisk after saving some space ### \n\n')
    overwrite_unneeded_fields(TraceList,rmFieldList);
    
    if o.truncateFeatureTS
        fprintf(1,'Shortening amax matrices ... ')
        ntr = numel(TraceList.eq.m);
        for itr = 1:ntr
            amax = TraceList.fb.amax{itr};
            pga  = TraceList.pga.tsi{itr};
            pgv  = TraceList.pgv.tsi{itr};
            pgd  = TraceList.pgd.tsi{itr}; % whos amax pga pgv pgd
            if ~isempty(amax); TraceList.fb.amax{itr} = amax(:,1:nsnpFeatures); end
            if ~isempty(pga);  TraceList.pga.tsi{itr} = pga (1:nsnpFeatures)  ; end
            if ~isempty(pgv);  TraceList.pgv.tsi{itr} = pgv (1:nsnpFeatures)  ; end
            if ~isempty(pgd);  TraceList.pgd.tsi{itr} = pgd (1:nsnpFeatures)  ; end
        end
        fprintf(1,' done.\n')
        TraceList.printObjectSize;
    end
    save(TraceListFullName ,'TraceList')
    save(SkipListFullName ,'SkipList')
    
    % Finish log-file
    unix(['echo "And run through it did." >> ' logFileName]);
    unix(['echo "Other scripts that were used in this run: " >> ' logFileName]);
    unix(['cat *.m >> ' logFileName]);
    diary off
    unix(['echo "########################## Matlab Output ###############################" >> ' logFileName]);
    unix(['cat ',diaryFileName,' >> ' logFileName]);
    unix(['rm ',diaryFileName]);
    unix(['bzip2 -f ',logFileName]);
    
    %unix('rm trList_tmp.mat');
    %unix(sprintf('rm %s',TraceListFullName_mod))
end