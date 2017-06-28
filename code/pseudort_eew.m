% Former fbank_paramInference.m script.
% Runs all records from traceList data set through a pseudo real-time EEW
% environment and simulates an EEW system. Makes various kinds of source
% estimates and ground motion predctions.
% 
% menandrin@gmail.com, 161219

clear all

% THINGS IMPLEMENTED IN fbank_paramInference_finder.m, but not here
% . test if amps in 72282711.BK.CVS.HNZ.00.sac are correct
% . ...
% . ...

configFileName = 'gba/gba_i39_TR3.m'; 
%configFileName = 'tpgx_i39.txt';
%configFileFullName = sprintf('~/programs/seismo/matlab/m_stats/config/%s',configFileName);
configFileFullName = sprintf('~/programs/seismo/matlab/trlist/config/%s',configFileName);

global fc o out snpLength snippet W mm rr minimax 
global fMode fOrder ftSize iN vp vs dataSetNames colours zDefault

addpath(genpath('gba/'))
addpath(genpath('m_tools/'))
addpath(genpath('../../trlist/'))
addpath(genpath('../../base/'))
addpath(genpath('../../../../matlab/'))

isMasterScript=true;    % Flag used when more than one inference script is 
                        % run in 'pathetically parallel' mode. Check "2 
                        % ways to run this script" section at beginning of script.


%% 2 ways to run this script:
%  i)  As a stand alone script in matlab,
%      for smallish projects and for computing reg. coeffs in parfor loop
%  ii) To speed things up, run multiple copies ('clones') simultaneously.
%      e.g. for parameter inference on large target-list.
%      - run ./make_fbpi_clone.sh   (deletes tmp/-directory & old clones and
%                                   copies this script to 
%                                   fbank_paramInference_clone.m, which has 
%                                   isMasterScript flag set to false)
%      - start screen session and run both fbank_paramInference.m and a
%        'reasonable' number of the fbpi_clone.m scripts.
%
%      8ung: do not run other copies of this script than the clones at the 
%            same time, because otherwise tmp-files will be overwritten.



%% NOTES

% PARAMETER INFERENCE CODE STRUCTURE
% for ieq = 1:neq                                         % Loop over all eqs
%   -->  Separate target and training data, prepare eq-wise GMP
%   
%     for it = 1:nt                                     % Loop over all time points
%       for itrg = 1:ntrig                                 % Loop over all triggered records at this time point
%       --> Estimate GMP & Src params
%       end
%     end
% end
%
% PEAK AMPLITUDE FIELDS
% amax                  fb matrix, on HP signals
% pa/pv/pd              time series, depending on options in fbank_wformProc_xx.m
%                       from BP or HP signals
% ppg[avd]/spg[avd]     peak amps 
% intensity
%
%  OUT-DIRECTORY NAMES
%  outDirFullName:      out-directory for a particular filter-bank run, e.g. 
%                       fbout/i32/causal_0p5     
%  projectDirFullName:  specific project within outDirFullName, using a
%                       specific set of parameters, e.g.
%                       fbout/i32/causal_0p5/SouthNapaEQ
% TASKS & IDEAS
%  - real time azimuth estimation by maximising low freqs on horizontals?
%
% GENERAL
%  For EEW GMPs: use HP filtered traces & 
%  For onset studies, Meier et al., 2016, GRL, use BP filtered traces
%  GbA uses: amax: HP filtered vel on individual comps
% 
%  ONSET TIMES: Only Japanese and Socal data have absolute origin times. 
%  For NGA data I can only estimate t0 from tppx and theor. travel times.
%
% Use add2blackList(1264,zList,'picked on foreshock','blackList_i29.mat')
% to add trace to blackList
%
% ISSUES
%  - dont use nga traces for estimating coeffs of 1st and 9th band
%  - Why so many events with only one record? check SkipList...  SNR?
%  - Jeremy has sent an email with comments on the variance issue for the weighted regression





%% %%%%%%%%%%%%%%%%%%
%  0. PREPARATIONS  %
%%%%%%%%%%%%%%%%%%%%%

% Load all settings and options from CONFIG file  . . . . . . . . . . . . . 
fprintf(1,['\n\n\tUsing config file: ',configFileName,'\n\n'])
fileId = fopen(configFileFullName);
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
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .


if isMasterScript && o.saveOut
    fprintf(1,'Note: if you intend to run clones, produce a clone file BEFORE you start this master script (with make_fbpi_clone.sh)\n')
    pause
end

tens         = linspace(1e1,1e5,1e4);
hundreds     = linspace(1e2,1e5,1e3);
thousands    = linspace(1e3,1e6,1e3);
tenthousands = linspace(1e4,1e6,1e2);

scriptName = [mfilename,'.m'];  % Script file name
pid        = feature('getpid'); % Process ID

allbands = (1:9)';
nsimbak  = nsim;
%nsimhalf = floor(nsim/2);  


% Parameter grid
mmin = minimax.mmin;
mmax = minimax.mmax;
mm   = mmin:.2:mmax;        % Grid
%rr   = 0:.1:2;              % log-r grid: log(1km): .. : log(100km)
rr   = log10([1,3:3:12,15:5:50,60:10:100]);

colours  = {'k','r','c','y','m','b',[0 .4 0],[255/256 165/256 0],[238/256 130/256 238/256],'k','m','c','g','y','b','r'};



%% %%%%%%%%%%%%%%
%  1. DATA SET  %
%%%%%%%%%%%%%%%%%

%  1.1 Load data sets
%  ==================
outDirName = sprintf('%s%s',fMode,outDirAppendix);
fprintf(1,['\n\nLoading filterBank.m-output from\t <basenames>/',outDirName,'/ ..  \n'])

outDirFullName       = sprintf('/scratch/memeier/fbout/i%i/%s/',iN,outDirName);
finderOutDirFullName = '/scratch/memeier/finder/out/';
projectDirFullName   = strcat([outDirFullName,out.projectName,'/']);
out.outDirFullName   = outDirFullName;

if ~exist(outDirFullName,'dir')     && o.saveOut; unix(['mkdir -p ',outDirFullName]);     end
if ~exist(projectDirFullName,'dir') && o.saveOut; unix(['mkdir -p ',projectDirFullName]); end

dataSetNames = strrep(dataSetBaseNames,'/out/',sprintf('/out/i%i/%s/',iN,outDirName));

TraceList = import_trLists(dataSetNames,'traceList');
fc        = TraceList.prop.fc;
snpLength = TraceList.prop.snpLength;
if o.useSkipList; SkipList = import_trLists(dataSetNames,'skipList'); end

nbands  = numel(fc(:,1));

fprintf(1,'Setting dist.flt=dist.hyp for traces that have no dist.flt-entry.\n')
idx=find(TraceList.dist.flt==0);
TraceList.dist.flt(idx) = TraceList.dist.hyp(idx);


tmpList = TraceList.selectSubList(2353:6000);



%% %%%%%%%%%%%%%%%%%%
% X. USE BLACKLIST  %
%%%%%%%%%%%%%%%%%%%%%
% Remove traces listed in blackList (starting i39: has already been done in wform proc script).
if o.useBlackList
	
    opts.process                = 1;
    opts.intMode                = 'afterPx';
    opts.scp_wforms             = true;
    opts.plot_blackListedTraces = 0;
    blackListFullName           = sprintf('~/programs/seismo/var/blackLists/blackList_i%i.mat',iN);

    TraceList = delete_blackListed_records(TraceList,blackListFullName,opts);
end


%% %%%%%%%%%%%%%%%%%%%
% X. OVERWRITE PICKS %
%%%%%%%%%%%%%%%%%%%%%%
% Uses picks stored in pxList and then recomputes Pa, Pv & Pd vectors as
% well as snr.
if o.usePxList
    
    % Set parameters for recomputing peak amps and Co.
    px              = load_picker_settings('production_short');
    opts.process    = 1;
    opts.intMode    = 'afterPx';
    opts.scp_wforms = true;
    opts.plotNewPx  = false;
    opts.ntap       = px.Param.ntap;
    opts.fLow_prefilt = px.Param.fLow_prefilt;
    opts.fOrder     = fOrder;
    opts.fMode      = fMode;
    pxListFullName  = sprintf('~/programs/seismo/var/pxLists/pxList_i%i.mat',iN);

    [TraceList] = overwritePx(TraceList,pxListFullName,px,opts);   
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X. RM OUTLIERS AND CLIPPERS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if o.excludeOutliers    
    isOut     = cell2mat(cellfun(@(x) sum(strcmp(x,'outlier')),TraceList.comment,'uniformOutput',0))>0;
    fprintf(1,sprintf('Removing %i/%i traces from z/e/nLists because they have been classified as OUTLIERS.\n',sum(isOut),numel(TraceList.eq.m)))
    %outList  = TraceList.selectSubList(find(isOut));
    TraceList = TraceList.selectSubList(find(~isOut));
end
if o.excludeClippers
    isClipped = cell2mat(cellfun(@(x) sum(strcmp(x,'clipped')),TraceList.comment,'uniformOutput',0))>0;
    fprintf(1,sprintf('Removing %i/%i traces from z/e/nLists because they have been classified as CLIPPED.\n',sum(isClipped),numel(TraceList.eq.m)))
    TraceList = TraceList.selectSubList(find(~isClipped));
end


 
% Time ....................................................................
tms     = snpList*snpLength;     % Multi station inference times
nt      = numel(tms);               % No. of time steps after t0
nsnpmax = max(snpList);
fprintf(1,'8ung: nsnpmax is used differently in different context. make sure it is uniquely defined. check. now.\n')
gmpOpts.tms = tms;

% Intensity  ..............................................................
mmiThresholdVect = TraceList.prop.immiThresholds; 
nit               = numel(mmiVectPex);  %OLD: nit= numel(mmiThresholdVect);
% .........................................................................


% Make sure amax entries are long enough for paramEstimation, i.e. at least
% as long as nsnpmax. SHOULD BECOME OBSOLETE IN FUTURE RUNS
tic
fprintf(1,'Make all amax- and pa/v/d entries have at least nsnpmax samples ...')
nsnpList = cell2mat(cellfun(@(x) size(x,2), TraceList.fb.amax,'uniformOutput',false));
if min(nsnpList)<nsnpmax &&~isempty(TraceList.fb.amax{1})
    
    ntr  = numel(TraceList.eq.m);
    for itr = 1:ntr
        
        print_iteration_numbers(itr,ntr,'tenthousands')
        
        ns = nsnpList(itr);
        if ns<nsnpmax; A                      = TraceList.fb.amax{itr};
                       TraceList.fb.amax{itr} = [A,repmat(A(:,end),1,nsnpmax-ns)];    
        end
        
        na = numel(TraceList.pga.tsi{itr});
        nv = numel(TraceList.pgv.tsi{itr});
        nd = numel(TraceList.pgd.tsi{itr});
        if na<nsnpmax; pa                     = TraceList.pga.tsi{itr};
                       TraceList.pga.tsi{itr} = [pa,repmat(pa(end),1,nsnpmax-na)];
        end
        if nv<nsnpmax; pv                     = TraceList.pgv.tsi{itr};
                       TraceList.pgv.tsi{itr} = [pv,repmat(pv(end),1,nsnpmax-nv)];
        end
        if nd<nsnpmax; pd                     = TraceList.pgd.tsi{itr};
                       TraceList.pgd.tsi{itr} = [pd,repmat(pd(end),1,nsnpmax-nd)];
        end
    end
end
toc
% Check
% nsnpList = cell2mat(cellfun(@(x) size(x,2), TraceList.fb.amax,'uniformOutput',false));
% na = min(cell2mat(cellfun(@(x) size(x,2), TraceList.pa,'uniformOutput',false)));
% nv = cell2mat(cellfun(@(x) size(x,2), TraceList.pv,'uniformOutput',false));
% nd = cell2mat(cellfun(@(x) size(x,2), TraceList.pd,'uniformOutput',false));
% if ~isequal(min(nsnpList),nsnpmax); error('sth went wrong. check.'); end


%  1.3 Split <traceList> into <zList> & <hList>
%  ============================================
fprintf(1,'\nSplitting TraceList into z-, e- & nList ... ')

idx_Z = find(strcmp(TraceList.station.ocode,'Z'));
idx_E = find(strcmp(TraceList.station.ocode,'E'));
idx_N = find(strcmp(TraceList.station.ocode,'N'));

zList = TraceList.selectSubList(idx_Z);
eList = TraceList.selectSubList(idx_E);
nList = TraceList.selectSubList(idx_N);
%clear TraceList
 
nz = numel(zList.eq.m);
ne = numel(eList.eq.m);
nn = numel(nList.eq.m);
if ( (nz~=ne) || (nz~=nn) )
    error('z-, e- & n-lists have different length. go find the bug. Now.\n')
end
fprintf(1,' done.\n')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%
% X. REMOVE DOUBLE ENTRIES  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'Removing doubly contained records ...') 
doubleIdxList = find_double_records(zList);
uniqueIdx     = setdiff(1:nz,doubleIdxList);
zList         = zList.selectSubList(uniqueIdx);
eList         = eList.selectSubList(uniqueIdx);
nList         = nList.selectSubList(uniqueIdx);
nz = numel(zList.eq.m);
ne = numel(eList.eq.m);
nn = numel(nList.eq.m);
fprintf(1,sprintf(' %i cases found. done.\n',numel(doubleIdxList)));




% Save some parameters in zList.prop
zList.prop.gmpOpts = gmpOpts;
zList.prop.tms     = tms;
zList.prop.configFileName = configFileName;


if o.useSkipList 
    idx_Z     = find(strcmp(SkipList.station.ocode,'Z'));
    zSkipList = SkipList.selectSubList(idx_Z);
end


% Clone nList to hList and modify amax values to arithmetic means of eList.fb.amax and nList.fb.amax.
fprintf(1,'Combining eList & nList to hList ... ')
hList = nList.selectSubList(1:nz);
%dar = zeros(nz,1); dvr = zeros(nz,1);

%zList.intensity.zen_ts = zeros(nz,nsnpmax);
for iz = 1:nz
    if ~isempty(eList.fb.amax{iz})
        % Combine filter bank peak amps
        aE                = eList.fb.amax{iz}(:,1:nsnpmax); 
        aN                = nList.fb.amax{iz}(:,1:nsnpmax);
        hList.fb.amax{iz} = 1/2*(aE+aN);
    end
end
fprintf(1,' done.\n')


% Clone flt-field from nList to zList for NGA traces
ngaIdx = find(strcmp(zList.dataSetName,'ngawest1'));
if ~isempty(ngaIdx); 
    zList.station.filter(ngaIdx) = nList.station.filter(ngaIdx);
    fprintf(1,'Note: Writing horizontal flt-values (from nList) into zList.flt\n')
end



% Exclude traces with snr<snrMin
if o.useSnrMin
    fprintf(1,'This section may be outdated. Check SNR definition.\n')
    lgcZ_lowSnr = logical(zList.noise.snr<snrMin);
    lgcH_lowSnr = logical(hList.noise.snr<snrMin);
    oneLow      = (lgcZ_lowSnr | lgcH_lowSnr);
    bothLow     = (lgcZ_lowSnr & lgcH_lowSnr);
    
    fprintf(1,sprintf('\nOf %i records, \n\t%i have one component with SNR<%3.1f, \n\t%i have both components with SNR<%3.1f, \n\t%i have vertical components with SNR<%3.1f\n', ...
        nz,sum(oneLow),snrMin,sum(bothLow),snrMin,sum(lgcZ_lowSnr),snrMin))
    
    fprintf(1,sprintf('Removing the %i records with vertical SNR<%i\n',sum(lgcZ_lowSnr),snrMin))
    zList = zList.selectSubList(find(~lgcZ_lowSnr));
    eList = eList.selectSubList(find(~lgcZ_lowSnr));
    nList = nList.selectSubList(find(~lgcZ_lowSnr));
    hList = hList.selectSubList(find(~lgcZ_lowSnr)); 
end
nz          = numel(zList.eq.m);

% Construct name for eqs-structure that is unique to the compiled zList
nzhalf          = round(nz/2);
eqsNameAppendix = sprintf('r%ikm_r%ikm_r%ikm',round(zList.dist.hyp([1,nzhalf,nz])));
eqsFullName     = sprintf('%seqs/eqs_ntr%i_%s.mat',outDirFullName,nz,eqsNameAppendix);

if o.checks; do_list_check(zList,eList,nList); end
fprintf(1,'done.\n')



%% %%%%%%%%%%%%%%%%%%%%%%%%%
% S. COMPILE STATION LIST  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'\n\nImporting station lists\n')

% Construct unique file name 
idxMiddle  = floor(nz/2);
middleDist = round(zList.dist.hyp(idxMiddle));
stationListFullName = sprintf('/scratch/memeier/fbout/i%i/stations/unifiedList_ntr%i_rmean%i.mat',iN,nz,middleDist);
% A. Assign unique event index to each trace: zList.eq.idx .................
if logical(exist(stationListFullName,'file'))
    
    fprintf(1,'loading existing unidfied stationlist-file...')
    load(stationListFullName)
    nst         = numel(stationList.lat);
    fprintf(1,'done.\n') 
else
    optVs.plotVsProfile  = false;
    optVs.printVsProfile = false;
    optVs.plotVsMap      = false;
    optVs.vs30_default   = vs30_default;
    stationList          = compile_unified_stationList(zList,optVs);
    save(stationListFullName,'stationList')
end



fprintf(1,'Assign station indices... ')
nz = numel(zList.eq.m); %zList = zList.selectSubList(find(zList.eq.m>=4));
%zList.station.vs30 = zeros(nz,1);
flg_notFound = false(nz,1);
hasDouble    = false(nz,1);
doubleStIdx  = cell(nz,1);
ct=0;
for iz = 1:nz
    
    print_iteration_numbers(iz,nz,'tenthousands')
    %zList.printSingleTraceSummary(iz)
    [stIdx,idxMultiples] = find_stationList_index(zList.station.lat(iz),zList.station.lon(iz),zList.station.name{iz},stationList);

    if ~isempty(idxMultiples)
        stationList.region(idxMultiples)
        stName = unique(stationList.name(idxMultiples));
        idx    = find(strcmp(zList.station.name,stName));
        zList.fullName(idx)
        hasDouble(iz) = true;
    end
    if ~isempty(stIdx);  zList.station.idx(iz)   = stIdx;
                         thisVs30                = stationList.vs30(stIdx);
                         if thisVs30<=0; thisVs30 = vs30_default; end
                         zList.station.vs30(iz)  = thisVs30;
    else                 zList.station.idx(iz)   = nan;
                         zList.station.vs30(iz)  = vs30_default;
                         flg_notFound(iz)        = true;
    end
end
fprintf(1,'done.\n')
fprintf(1,sprintf('%i/%i traces had no list entry in stationList.\n',sum(flg_notFound),nz))
% idxDouble = find(hasDouble);
% nd        = numel(idxDouble);
% for id = 1:nd
%     zList.printSingleTraceSummary(idxDouble(id))
%     idx = doubleStIdx{idxDouble(id)};
%     zList.fullName(idx)
% end




%% %%%%%%%%%%%%%%%%%%%%%%%%
% 3. COMPILE EVENT LISTS  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if o.make_eqs_list
    
    % zList.eq.idx      <nz-by-1 >      event-index of each trace
    % eqs.eventId       <neq-by-1>      list of unique eqIds
    % eqs.traceId       <neq-by-1>      cell with traceList-indices of all
    %                                   traces of respective eqs.eventId,
    %                                   sorted wrt/ hypocentral distance
    % eqs.tpx           <neq-by-1>      cell with theoretical p-wave arrival
    %                                   times of traces in eqs.traceId
    
    fprintf(1,'\nCompiling list with all eqs ...\n')
    
    % A. Assign unique event index to each trace: zList.eq.idx .................
    if logical(exist(eqsFullName,'file'))
        
        fprintf(1,'loading existing eqIdList-file...')
        load(eqsFullName)
        zList.eq.idx = eqs.eqIdx;
        eqs          = rmfield(eqs,'eqIdx');
        fprintf(1,'done.\n')
       
    else
        
        if o.useSkipList; eqs = compile_eqs_structure(zList,zSkipList);
        else              eqs = compile_eqs_structure(zList,[]);
        end
        fprintf(1,'saving ... ')
        save(eqsFullName,'eqs')
        fprintf(1,'done.\n')
    end
    if exist('zSkipList','var'); clear zSkipList; end
    if exist('SkipList' ,'var'); clear SkipList;  end
    fprintf(1,'done\n')
    %if o.checks; check_if_eventIds_increase(eqs); do_eqsName_check(eqs,zList); end
    
    zList.prop.eqs = eqs;
    
    % Merge double entries from NGA & Socal ...............................
    fprintf(1,'Merging double entries from NGA & Socal. \n')
    eqs = merge_double_eqs_entries(eqs,'3031111','LANDERS_1992');
    eqs = merge_double_eqs_entries(eqs,'9108652','HECTOR_1999' );
    eqs = merge_double_eqs_entries(eqs,'3144585','NORTHR_1994' );
    neq = numel(eqs.eventId);
end

fprintf(1,'\n8UNG: eqs NOW CONTAINS zList-INDICES FOR ALL EARTHQUAKES. IF YOU CHANGE zList YOU NEED TO RECOMPUTE eqs\n')
isJapa = strcmp(zList.dataSetName,'kNet') |strcmp(zList.dataSetName,'kikNet');
isCali = strcmp(zList.dataSetName,'scsn') |strcmp(zList.dataSetName,'scsnPx');
isNgaw = strcmp(zList.dataSetName,'ngawest1');
isWnch = strcmp(zList.dataSetName,'wenchuan');
%print_eqsList(eqs,'/scratch/memeier/data/composite/eqsLists/i38/new/eqsList_Mge4.txt')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  --> LOTS OF OUTDATED STUFF MOVED TO APPENDIX APPENDIX APPENDIX    <--  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zList.prop.gmThresh.mmi.pex = mmiVectPex;
zList.prop.gmThresh.mmi.tex = zList.prop.immiThresholds;
% FUTURE: zList.prop.gmThresh.mmi.tex --> already populated. CHECK.



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE TRACELIST INTO NETCDF DATABASE % (section copied to pseudort_eew.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if o.tracelist2dataBase
    
    slashIdx    = regexp(configFileName,'/');
    nc.fileName = strrep(configFileName(slashIdx+1:end),'.m','.nc');
    nc.fullName  = sprintf('%s%s',nc.path,nc.fileName);
    
    tracelist2netcdf_gbaTraining(zList,hList,nc.fullName,nc.nsnippet)
    print_asciiTable_simple(zList,strrep(nc.fullName,'.nc','.txt'));

    configCopyFullName = strrep(nc.fullName,'.nc','.config');
    ucmd = sprintf('cp -p %s %s',configFileFullName,configCopyFullName);
    unix(ucmd);
    %save(strrep(ncFullName,'.nc','_trList.mat'),'zList')

    % Read Netcdf file (as test)
    %vardata    = ncread(ncFullName,'z');
    %zFullNames = ncread(ncFullName,'zFullName');
    %nFullNames = ncread(ncFullName,'h1FullName');
    %eqIdx      = ncread(ncFullName,'eqIdx');
end



%% %%%%%%%%%%%%%%%%%%%%%%%%
% 5. PARAMETER INFERENCE  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5.1.2 Prepare parallelisation
procEqListFullName = strcat(['~/programs/seismo/matlab/projects/eew/tmp/',num2str(neq),'eq_',num2str(nz),'traces.txt']);

if isMasterScript

	% Check if there are older .mat files in the project directory 
    if o.saveOut
        fNames    = dir(projectDirFullName);
        fNames    = {fNames.name}';
        nMatFiles = sum(cellfun(@(x) ~isempty(x), regexp(fNames,'*.mat','match')));
        if nMatFiles~=0
            fprintf(1,'8ung, there are existing *.mat-files in the project-out-directory. Delete or move them manually before continuing.\n')
            pause
        end
    end

    % Start a list containing indices of eqs (ieq) that have already been processed
    if exist('tmp','dir'); unix('rm -rf tmp/'); end
    unix('mkdir tmp');
    unix(['echo ''0 0'' > ',procEqListFullName]); % Create blank file
    
else
    % Wait until master script has initiated the allTarget-list
    flg_ready = false;
    while ~flg_ready
        if ~exist(procEqListFullName,'file')
            fprintf(1,'procEqList not ready yet. Waiting for masterScript to start procEqList. Trying again in 1min.\n')
            pause(60)
        else
            flg_ready=true;
        end
    end
end





%% 5.2 Loop over all earthquakes
% ==============================            check LOOP_OVERVIEW_PSEUDOFCT.m

% Make source and ground motion estimations/predictions in pseudo-real-time
% for all earthquakes. Save estimations based on single waveforms in 
eqs.rte   = cell(neq,1);
ct_eqProc = 0;
ct_eqs    = 0;

if o.testIndexAssignment; test_index_assignment(zList,eqs,tms,rmaxTarget,snpLength); end

%dsList        = zList.dataSetName(cellfun(@(x) x(1), eqs.traceId));
%is_wch_or_nga = logical(strcmp(dsList,'ngawest1') |strcmp(dsList,'wenchuan'));
%neq_targets   = sum(eqs.m>=mminTarget &eqs.m<mmaxTarget &~is_wch_or_nga);
%neq_targets   = sum(eqs.m>=mminTarget &eqs.m<mmaxTarget);
%neq_targets   = 1;

% Define target list
isTargetEq  = strcmp(eqs.name,'15481673'); % LaHabra
%isTargetEq  = strcmp(eqs.name,'72282711'); % South Napa
neq_targets = sum(isTargetEq);  


%profile on
for ieq = 1:neq

    % Check if eq has already been processed
    %ucmd         = strcat(['grep ''^',num2str(ieq),'\ '' ', procEqListFullName]);
    %[status,hit] = unix(ucmd);
    hit = [];
    
    % If not, and if the mangitude is big enough
    %meetsEqTargetCriterion = eqs.m(ieq)>=mminTarget &eqs.m(ieq)<mmaxTarget &~is_wch_or_nga(ieq);
    %meetsEqTargetCriterion = eqs.m(ieq)>=mminTarget && eqs.m(ieq)<mmaxTarget;
    %meetsEqTargetCriterion = strcmp(eqs.name{ieq},'15481673');
    %if meetsEqTargetCriterion    
    
    if isTargetEq(ieq)

        %if (isempty(hit) && eqs.m(ieq)>=mminTarget && eqs.m(ieq)<mmaxTarget)
        %if  (eqs.m(ieq)>5.169 &eqs.m(ieq)<5.171)                  % Borrego Springs 2016
        ct_eqProc = ct_eqProc+1;
         
        % Add eq to procEqList and save it
        %printLine = sprintf('%d \t %d',ieq,pid);
        %ucmd      = strcat(['echo ''',printLine,''' >> ',procEqListFullName]);
        %unix(ucmd);
        
        %nsnpmax = numel(tms);
        %nt    = numel(tms);

        fprintf(1,sprintf('\n\n\nEARTHQUAKE %s, M%3.1f on %s\n',eqs.name{ieq},eqs.m(ieq),eqs.date{ieq}))
        % 5.2.1 Read eq-info, separate target & training data, ...
        % ========================================================
        
        % INDEX NAMING CONVENTION .........................................................................
        % Two levels of traceLists, each of which has indices that point to entries in list
        % 1. zList/hList    . idxTargets,        all records of the same earthquake
        %                   . idxTraining,       all records NOT of the same earthquake
        %                   . idxTriggers,       all records that at it^th time step have triggered
        %                   . idxCurrentTrigger, record for which rte is being made at current time step
        % 2. zTargetList    . iiTargets,         same as idxTargets but relative to zTargetList
        %                   . iiCurrentTrigger,  same as idxCurrentTrigger but relative to zTargetList
        %                   . iiTriggers,        same as idxTriggers but relative to zTargetList
        %
        % If you want to find the zList-index of the itarget^th entry in zTargetList, use idx = idxTargets(itarget);
        % .................................................................................................
        
        
        
        %%%%%%%%%%%% SEPARATE TARGET & TRAINING DATA  %%%%%%%%%%%%%%%%%%%%%
        idxTargets  = eqs.traceId{ieq};
        idxTraining = setdiff(1:nz,idxTargets)'; % traceList index of all traces not from same event
        
        
        %%%%%%%%%%%%%%%%%% PREPARE TARGET DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(1,'Preparing target data ...'); tic
        isValidTarget = zList.dist.hyp(idxTargets)<rmaxTarget;
        idxTargets    = idxTargets(isValidTarget);
        zTargetList   = zList.selectSubList(idxTargets);
        hTargetList   = hList.selectSubList(idxTargets);
        rTargets      = zTargetList.dist.hyp;
        ntargets      = numel(idxTargets); %ntargets   = numel(rTargets<=rmaxTarget);
        
        tppx       = zTargetList.px.p.t;
        t_firstPpx = min(tppx);
        
        if o.writeGbaAsciiOut
            gbaAscii.m.mu           = zeros(ntmax_ascii,ntargets);
            gbaAscii.m.var          = zeros(ntmax_ascii,ntargets);
            gbaAscii.r.mu           = zeros(ntmax_ascii,ntargets);
            gbaAscii.r.var          = zeros(ntmax_ascii,ntargets);
            gbaAscii.pgaPred.mu     = zeros(ntmax_ascii,ntargets);
            gbaAscii.pgaPred.var    = zeros(ntmax_ascii,ntargets);
            gbaAscii.pgvPred.mu     = zeros(ntmax_ascii,ntargets);
            gbaAscii.pgvPred.var    = zeros(ntmax_ascii,ntargets);
            gbaAscii.az             = zeros(ntmax_ascii,ntargets,nbands);
            gbaAscii.ah             = zeros(ntmax_ascii,ntargets,nbands);
            gbaAscii.tppx           = zeros(ntargets,1);
        end
        
        
        % Extract target lists, initiate ssrte structures, make sure
        % amax-matrices are long enough.
        % Extracting individual targets from main traceList (zList) is the 
        % most expensive step. It saves a lot of time to extract them for 
        % each site of an earthquake to an array once in the beginning, 
        % rather than extracting it at each time step for each triggered site
        fprintf(1,'.. Extract the targetList for each target site and initiate ssrte-structures... ')
        zTargetArray     = cell(ntargets,1);
        hTargetArray     = cell(ntargets,1);
        ssrteArray       = cell(ntargets,1);
        ssrte            = initiate_ssrte(rteMethodList,nt,nit,gmpOpts);
        ssrte.t.firstPpx = t_firstPpx;
        ssrte.dy         = cell(nt,1);
        for itarget = 1:ntargets
            print_iteration_numbers(itarget,ntargets,'hundreds')
            ssrteArray  {itarget} = ssrte;
            zTargetArray{itarget} = zTargetList.selectSubList(itarget);
            hTargetArray{itarget} = hTargetList.selectSubList(itarget);
        end
        fprintf(1,'done\n.')

        if (numel(unique(zTargetList.eq.m))>1); fprintf(1,'\t8UNG: more than one magnitude in target-list ... sth wrong with event-ID\n'); end
        
        % Eq meta info
        if ~isempty(idxTargets)
            eqLat  = zTargetList.eq.lat (1);
            eqLon  = zTargetList.eq.lon (1);
            eqZ    = zTargetList.eq.z   (1);
            eqDate = zTargetList.eq.date{1};
            m      = zTargetList.eq.m   (1);
            eqName = zTargetList.eq.name{1};
        end
        
        targetIsJapa = strcmp(zTargetList.dataSetName{1},'kNet') |strcmp(zTargetList.dataSetName{1},'kikNet');
        if targetIsJapa; ngaParams.region = 3;
        else             ngaParams.region = 0;
        end
        
        
        %%%%%%%%%%%%%%%%%% PREPARE TRAINING DATA  %%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(1,'Preparing training data ...'); tic

        % All training criteria must be relative to zList
        %isValidTraining = logical(zList.dist.hyp(idxTraining)<=rmaxTraining &zList.noise.snr(idxTraining)>=snrMin &~isNgaw(idxTraining) &~isWnch(idxTraining));
        isValidTraining = logical(zList.dist.hyp(idxTraining)<=rmaxTraining &zList.noise.snr(idxTraining)>=snrMin);
        idxTraining     = idxTraining(isValidTraining);
        zTrainingList   = zList.selectSubList(idxTraining);
        hTrainingList   = hList.selectSubList(idxTraining);
        ntraining       = numel(zTrainingList.eq.m);
        [noproblem]     = check_training_set(zTrainingList,hTrainingList);
        % CHECK: sum(strcmp(zTrainingList.comment,'clipped'));
        
        
        
        % TRAIN-STRUCTURE: The training data for the GbA will be passed on to estimate_params_GbA.m 
        % as <train> structure containing observed peak filer bank amplitudes, magnitudes, hypocentral distances, 
        % peak GM values. Some computational time can be saved by only
        % passing on the peak amplitudes for 
        % Read amp-matrices for each snippet
        % Precompute snippet indices that will be processed (saves time) ...
        fprintf(1,'. Identify snippets that will be used ...')
        tic
        usedSnippets = [];
        for it = 1:nt
            
            % Find triggered stations
            tnow    = tms(it);                  % Time since origin time
            tavail  = tnow - tppx;              % How much time since P-onset on each target record?
            hasTrig = logical(tavail>=0);
            %hasTrig = logical(tavail>=0 &rTargets<=rmaxTarget);
            ntrig    = sum(hasTrig);
            if ntrig>0
                nsnpavail = round(tavail(hasTrig)/snpLength);
                nsnpavail(nsnpavail==0) = 1;
                usedSnippets = [usedSnippets; nsnpavail];
            end
        end
        usedSnippets = sort(unique(usedSnippets));
        fprintf(1,[' ',num2str(nsnpmax-numel(usedSnippets)),'/',num2str(nsnpmax),' snippets not used ... '])
        fprintf(1,'done. '); toc
        
        
        train.m   = zTrainingList.eq.m;
        train.r   = zTrainingList.dist.hyp;
        train.pga = zTrainingList.pga.pns.ampzen;
        train.pgv = zTrainingList.pgv.pns.ampzen;
        train.mmi = zTrainingList.mmi.pns.ampzen;
        
        % Now save all snippets that will be used into zamps and hamps, which will go into
        % the <train> structure which contains the training data used by GbA
        fprintf(1,'. Extracting amplitudes and saving them in training-structure ...')
        zamps = cell(1,nsnpmax);
        hamps = cell(1,nsnpmax);
        for isnp = 1:nsnpmax
            print_iteration_numbers(isnp,nsnpmax,'tens')
            if ismember(isnp,usedSnippets)
                zamps{isnp} = cell2mat(cellfun(@(x) x(:,isnp)', zTrainingList.fb.amax,'UniformOutput',0));
                hamps{isnp} = cell2mat(cellfun(@(x) x(:,isnp)', hTrainingList.fb.amax,'UniformOutput',0));
            end
        end
        fprintf(1,' done. '); toc
        
        
        
        %%%%%%%%%%%%%%%%%%%% MISCELLANEOUS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Write noise amplitudes into finder-files for all sites. Sites
        % that trigger during the course of an event will be overwritten
        %         if o.writeFinderFiles
        %             finderMat0 = zeros(ntargets,3);
        %             for itarget = 1:ntargets
        %                 logPaNoise            = log10(100*zTargetList.noise.acc{itarget}(3));
        %                 finderMat0(itarget,1) = zTargetList.station.lat(itarget);
        %                 finderMat0(itarget,2) = zTargetList.station.lon(itarget);
        %                 finderMat0(itarget,3) = logPaNoise;
        %             end
        %             
        %             finderObsFullName  = sprintf('%snew/%s/paObs/' ,finderOutDirFullName,eqName);
        %             finderPredFullName = sprintf('%snew/%s/paPred/',finderOutDirFullName,eqName);
        %             if ~exist(finderObsFullName ,'dir'); unix(sprintf('mkdir -p %s',finderObsFullName)) ; end
        %             if ~exist(finderPredFullName,'dir'); unix(sprintf('mkdir -p %s',finderPredFullName)); end
        %         end
        
        % Plot amplitude normalised traces of this event
        if o.plot_all_traces_of_eq; ht = plot_triggered_wforms(idxTargets,zList,tms,opts); end
        
        if o.useWt; wti  = Wti(idxTraining);
                    wtn  = Wtn(idxTraining);     % used for nonlinreg
                    wti0 = ones(size(idxTraining,1),1);
        else        wti0 = ones(size(idxTraining,1),1);
        end
        
        if ~o.quiet
            fprintf(1,['\n............................................................................................ \n'])
            fprintf(1,sprintf('New event, nr. %i/%i; m%3.1f %s;\t%i vert. traces from this event, %i with R<= %3.1fkm: \n', ...
                ct_eqProc,neq_targets,eqs.m(ieq),eqs.date{ieq},ntargets,numel(rTargets<=rmaxTarget),rmaxTarget))
            for itarget = 1:ntargets; fprintf(1,'  %s %5.1fkm\n',zTargetList.fullName{itarget},zTargetList.dist.hyp(itarget)); end
            
            if o.useSkipList;
                Rskip         = eqs.Rskip{ieq}; 
                nskip         = numel(Rskip);
                RskipProximal = Rskip(Rskip<rmaxTarget);
                nskipProximal = numel(RskipProximal);
                fprintf(1,sprintf('  %i skipped traces, %i with R<= %3.0fkm;',nskip,nskipProximal,rmaxTarget))
                fprintf(1,sprintf('%3.0f, ', Rskip))
                fprintf(1,sprintf('km\n'))
            end
        end

        % Initiate rte sturcture (used to be for likelihood = f(t))
        %rte                  = initiate_rte(rteMethodList,ntargets,nt,nit,gmpOpts);
        rte.idxTargets       = idxTargets;
        rte.stNames          = zTargetList.station.name;     % --> for back-checks
        rte.ntrig            = zeros(1,nt,'single');
        rte.nsecperstat      = zeros(ntargets,nt,'single');
        rte.ka13.pd.val      = zeros(ntargets,nt,'single');
        rte.ka13.m.hat       = zeros(1,nt,'single');
        rte.ka13.m.std       = zeros(1,nt,'single');
        rte.loc.rand.dr_hyp  = zeros(ntargets,nt,'single');
        rte.loc.rand.dr_epi  = zeros(ntargets,nt,'single');
        rte.loc.rand.std     = zeros(1,nt,'single');
        rte.loc.rand.lat.hat = zeros(1,nt,'single');
        rte.loc.rand.lon.hat = zeros(1,nt,'single');
        rte.gba.m.hat        = zeros(1,nt,'single');
        rte.gba.m.hat2       = zeros(1,nt,'single');
        rte.gba.m.std        = zeros(1,nt,'single');
        rte.gba.m.std2       = zeros(1,nt,'single');
        rte.t.firstPpx       = t_firstPpx;
        
        lmmi  = cell(ntargets,1);
        nsim  = nsimbak;
        bands = allbands;
        mHat  = 0;            % MLE w/o using r0, but changing <nsim> and <bands> when M_MLE>5.5
        
        %rte.empGM_ss.opts.plot=1;
        if sum(gmpOpts.plotHat)>0; gmpOpts.plotGMP=1;
        else                       gmpOpts.plotGMP=0;
        end
        
         
        % Loop over all time steps ..............................................................
        for it = 1:nt
        
            % Find triggered stations
            tnow          = tms(it);                  % Time since origin time
            tavail        = tnow - tppx;              % How much time since P-onset on each target record?
            hasTrig       = logical(tavail>=0);
            %hasTrig       = logical(tavail>=0 &rTargets<=rmaxTarget);
            ntrig         = sum(hasTrig);
            rte.ntrig(it) = ntrig;
            
            fprintf(['\nAt t=',num2str(tnow,'%4.1f'),'sec, ',num2str(ntrig),' triggers  '])
            
            % Limit max no. of triggers to ntrigMax in the interest of speed
            if ntrig>ntrigMax; ntrig=ntrigMax; end
            
            % Make list of idxTargets for the ntrig records with the most data available
            [~,iiTriggers] = sort(tavail,'descend');
            iiTriggers     = iiTriggers(1:ntrig);
            idxTriggers    = idxTargets(iiTriggers);
            % tmp1=zList.px.p.t(idxTriggers); tmp2=sort(zTargetList.px.p.t); isequal(tmp1,tmp2(1:ntrig))
            
            % if o.writeFinderFiles; finderObsMat=finderMat0; finderPredMat=finderMat0; end
            
            %% Onsite/single-station estimates on all currently triggered stations
            % Loop over all stations that triggered at this point in time  .  .  .  .  .  .  .  .  .  .  .
            for itrig = 1:ntrig
                
                if ntrig<10; fprintf([num2str(itrig),' .. '])
                else        
                    if ismember(itrig,10:10:1e3); fprintf([num2str(itrig),' .. ']); 
                    end
                end
                idxCurrentTrigger = idxTriggers(itrig);
                iiCurrentTrigger  = iiTriggers(itrig);
                tavail1           = tnow - tppx(iiCurrentTrigger);     % available data on this record
                nsnpavail         = round(tavail1/snpLength);
                if nsnpavail==0; nsnpavail=1; end
                
                % Check index assignment for current trace
                isSameTrace = strcmp(zList.fullName{idxCurrentTrigger},zTargetList.fullName{iiCurrentTrigger});
                if ~isSameTrace; fprintf(1,'Something wrong with index assignment. Check.\n'); pause; end
                tavail2 = tnow - zList.px.p.t(idxCurrentTrigger);     % available data on this record
                if ~isequal(tavail1,tavail2); fprintf(1,'checkmenow.\n'); pause; end
                % idxCurrentTrigger = idxTriggers(itrig);                       % zList index of single triggered target trace
                % iiCurrentTrigger   = find(idxCurrentTrigger==idxTargets);       % Index of current target relative to zTargetList
                
                if o.verbose; fprintf(1,'\n Select target...'); end

                zTarget = zTargetArray{iiCurrentTrigger};
                hTarget = hTargetArray{iiCurrentTrigger};
                ssrte   = ssrteArray  {iiCurrentTrigger};
                ssrte.t.avail(it) = tavail1;
                
                vs30target = zTarget.station.vs30;
                %                 if vs30target>vs30_threshold; pgam = pgam_rock;
                %                                               pgvm = pgvm_rock;
                %                 else                          pgam = pgam_soil;
                %                                               pgvm = pgvm_soil;
                %                 end
                
                gmpOpts.figFullName = sprintf('~/programs/seismo/fig/i37/gmp/realtime/new/gmpp_eq%i_tr%i',ieq,idxCurrentTrigger);
                %gmpOpts.tms         = tms;
                
                % Assign r-estimate from independent distance constraint
                r_ctlg      = zTarget.dist.hyp(1);
                r0_tmp.mean = -1;
                if ntrig<3; r0_tmp.sigma = r0.sigma;
                else        r0_tmp.sigma = r0.sigma2;
                end
                
                while r0_tmp.mean<0; r0_tmp.mean = r0_tmp.sigma*randn(1,1)+r_ctlg; end
                
                train.az = zamps{nsnpavail};
                train.ah = hamps{nsnpavail};
                nzilt_z  = numel(find(train.az(:)==0));
                nzilt_h  = numel(find(train.ah(:)==0));
                
                if o.verbose; fprintf(1,'\n Estimating params (GbA)...'); end
                gbaMap = estimate_params_GbA(zTarget,hTarget,train,wti0,nsnpavail,bands,{'Z','H'},nsim,[],[]);
                
                if o.writeGbaAsciiOut
                    
                    if nsnpavail<=ntmax_ascii
                        gbaAscii.m.mu (nsnpavail,iiCurrentTrigger) = gbaMap.m.hat2.mean;
                        gbaAscii.m.var(nsnpavail,iiCurrentTrigger) = gbaMap.m.hat2.std^2;
                        gbaAscii.r.mu (nsnpavail,iiCurrentTrigger) = gbaMap.r.hat2.mean;
                        gbaAscii.r.var(nsnpavail,iiCurrentTrigger) = gbaMap.r.hat2.std^2;
                        
                        gbaAscii.pgaPred.mu (nsnpavail,iiCurrentTrigger) = gbaMap.pga.mean;
                        gbaAscii.pgaPred.var(nsnpavail,iiCurrentTrigger) = gbaMap.pga.var;
                        gbaAscii.pgvPred.mu (nsnpavail,iiCurrentTrigger) = gbaMap.pgv.mean;
                        gbaAscii.pgvPred.var(nsnpavail,iiCurrentTrigger) = gbaMap.pgv.var;
                        
                        gbaAscii.az(nsnpavail,iiCurrentTrigger,:) = zTarget.fb.amax{1}(bands,nsnpavail);
                        gbaAscii.ah(nsnpavail,iiCurrentTrigger,:) = hTarget.fb.amax{1}(bands,nsnpavail);
                        
                        gbaAscii.tppx(iiCurrentTrigger) = zTarget.px.p.t;
                    end
                end
                
                % SRC PARAM ESTIMATION   ..............................
                % Magnitude estimations
                ssrte.gba_ss.m.hat (it) = gbaMap.m.mpe;
                ssrte.gba_ss.m.hat (it) = gbaMap.m.mpe;
                ssrte.gba_ss.m.hat2(it) = gbaMap.m.hat2.mean;
                ssrte.gba_ss.m.llh {it} = gbaMap.m.lmn';
                lmmi{iiCurrentTrigger}  = gbaMap.m.lmn';
                rte.gba_ss.meta.mm      = gbaMap.mm;
                
                % Distance estimations
                ssrte.gba_ss.rh.hat (it) = 10^gbaMap.r.mpe;
                ssrte.gba_ss.rh.hat2(it) = gbaMap.r.hat2.mean;
                ssrte.gba_ss.rh.llh {it} = gbaMap.r.lmn;
                
                ssrte.dy{it} =gbaMap.dy;
                % GM OBSERVATION ......................................
                %immiObs = zTarget.intensity{1};
                % .....................................................
                
                
                % GM PREDICTION .......................................
                % 1. Predict IMMI directly
                % 2. Predict mHat/rHat; predict PGx; predict IMMI
                % 3. Take    mObs/rObs; predict PGx; predict IMMI
                % .....................................................
                
                % 1. Predict IMMI directly    .  .  .  .  .  .  .  .  .
                if ismember('empGM_ss',rteMethodList)
                    ssrte.empGM_ss.mmi.hat (it)  = gbaMap.mmi.mean;
                    ssrte.empGM_ss.mmi.var (it)  = gbaMap.mmi.var;
                    ssrte.empGM_ss.mmi.prct{it}  = gbaMap.mmi.prct;
                    ssrte.empGM_ss.pga.hat (it)  = gbaMap.pga.mean;
                    ssrte.empGM_ss.pga.var (it)  = gbaMap.pga.var;
                    ssrte.empGM_ss.pga.prct{it}  = gbaMap.pga.prct;
                    ssrte.empGM_ss.pgv.hat (it)  = gbaMap.pgv.mean;
                    ssrte.empGM_ss.pgv.var (it)  = gbaMap.pgv.var;
                    ssrte.empGM_ss.pgv.prct{it}  = gbaMap.pgv.prct;
                    %ssrte.empGM_ss.mmi.pEx(:,it) = get_propOfExceedence_ecdf(gbaMap.mmi.val,mmiThresholdVect);
                    ssrte.empGM_ss.mmi.pEx(:,it) = get_propOfExceedence_ecdf(gbaMap.mmi.val,mmiVectPex);
                end
                
                % 2. Use GbA mHat/rHat to predict PGx to predict IMMI
                mHat  = gbaMap.m.mpe;
                %rHat1 = gbaMap.r.mpe;       % mpe of log10(r), from maximising likelihood function
                %rHat2 = gbaMap.r.hat2.mean; % sample mean of nsim most similar traces's linear distances
                
                % NOTE: I STILL NEED TO CHECK IF MHAT2 IS SAME AS MHAT  --> isnot!
                mIn.mu  = gbaMap.m.hat2.mean;
                mIn.var = gbaMap.m.hat2.std^2;
                rIn.mu  = 10^gbaMap.r.hat2.mean;
                rIn.var = 999; % gbaMap.r.hat2.std^2;
                
                if ismember('gba_ss',rteMethodList)
                    fprintf('rIn.var is computed using sample variance of log-normal distribution. Make sure this is correctly translated into location uncertainty for GMPE... not yet done. do it. now.\n')
                    pause
                    if o.verbose; fprintf(1,'\n pGMP_MC for GbA...'); end
                    %nrnd_pgx=1e5; da=2e-2;
                    [pPGA,pgaEdges,pgaVect,pgaHat] = pGMP_MC(pgam,mIn,rIn,da,nrnd_pgx);
                    [pPGV,pgvEdges,pgvVect,pgvHat] = pGMP_MC(pgvm,mIn,rIn,da,nrnd_pgx);
                    %pgaCenters=pgaEdges(1:end-1)+da/2;  clf; plot(pgaCenters,pPGA); hold on; plot(pgaHat.val,pgaHat.p,'or')
                    
                    % Save PGx distributions
                    ssrte.gba_ss.pga.prct{it} = prctile(pgaVect,[1,5,16,50,84,95,99]);
                    ssrte.gba_ss.pga.hat (it) = pgaHat.val;
                    ssrte.gba_ss.pga.mean(it) = mean(pgaVect);
                    ssrte.gba_ss.pga.var (it) = var(pgaVect);
                    ssrte.gba_ss.pgv.prct{it} = prctile(pgvVect,[1,5,16,50,84,95,99]);
                    ssrte.gba_ss.pgv.hat (it) = pgvHat.val;
                    ssrte.gba_ss.pgv.mean(it) = mean(pgvVect);
                    ssrte.gba_ss.pgv.var (it) = var(pgvVect);
                    
                    % Turn PGx distributions into iMMI distributions
                    iMMI                       = pgx2immi(10.^pgaVect,10.^pgvVect);
                    pimmi                      = sampleVect2pdf(iMMI.immi,dimmi);
                    ssrte.gba_ss.mmi.hat (it)  = pimmi.mpe;
                    ssrte.gba_ss.mmi.mean(it)  = mean(iMMI.immi);
                    ssrte.gba_ss.mmi.var (it)  = var (iMMI.immi);
                    ssrte.gba_ss.mmi.prct{it}  = prctile(iMMI.immi,[1,5,16,50,84,95,99]);
                    ssrte.gba_ss.mmi.pEx(:,it) = get_propOfExceedence_ecdf(iMMI.immi,mmiVectPex);
                end

                
                if ismember('gban_ss',rteMethodList)
                    fprintf('rIn.var is computed using sample variance of log-normal distribution. Make sure this is correctly translated into location uncertainty for GMPE... not yet done. do it. now.\n')
                    pause
                    if o.verbose; fprintf(1,'\n pGMP_MC GBA_ss with NGAW2 GMPE...'); end
                    [pPGA,pgaEdges,pgaVect,pgaHat] = pGMP_MC_ngaw2('PGA',mIn,rIn,da,nrnd_pgx,vs30target,ngaParams);
                    [pPGV,pgvEdges,pgvVect,pgvHat] = pGMP_MC_ngaw2('PGV',mIn,rIn,da,nrnd_pgx,vs30target,ngaParams);
                    
                    % Save PGx distributions
                    ssrte.gban_ss.pga.prct{it} = prctile(pgaVect,[1,5,16,50,84,95,99]);
                    ssrte.gban_ss.pga.hat (it) = pgaHat.val;
                    ssrte.gban_ss.pga.mean(it) = mean(pgaVect);
                    ssrte.gban_ss.pga.var (it) = var (pgaVect);
                    ssrte.gban_ss.pgv.prct{it} = prctile(pgvVect,[1,5,16,50,84,95,99]);
                    ssrte.gban_ss.pgv.hat (it) = pgvHat.val;
                    ssrte.gban_ss.pgv.mean(it) = mean(pgvVect);
                    ssrte.gban_ss.pgv.var (it) = var (pgvVect);
                    
                    % Turn PGx distributions into iMMI distributions
                    iMMI                         = pgx2immi(10.^pgaVect,10.^pgvVect);
                    pimmi                        = sampleVect2pdf(iMMI.immi,dimmi);
                    ssrte.gban_ss.mmi.hat (it)  = pimmi.mpe;
                    ssrte.gban_ss.mmi.mean(it)  = mean(iMMI.immi);
                    ssrte.gban_ss.mmi.var (it)  = var (iMMI.immi);
                    ssrte.gban_ss.mmi.prct{it}  = prctile(iMMI.immi,[1,5,16,50,84,95,99]);
                    ssrte.gban_ss.mmi.pEx(:,it) = get_propOfExceedence_ecdf(iMMI.immi,mmiVectPex);
                end
                
                
                % 3. "PYTHIA": Use catalog values mCtlg/rCtlg to predict PGx to predict MMI
                %    . using CH07  (pythia_ss)
                %    . using NGAW2 (pythia2_ss)
                mIn.mu  = zTarget.eq.m;
                mIn.var = 0;
                rIn.mu  = zTarget.dist.flt;
                rIn.var = 0;
                
                if ismember('pythia_ss',rteMethodList)
                    if o.verbose; fprintf(1,'\n pGMP_MC for Pythia...'); end
                    %nrnd_pgx=1e5; da=2e-2;
                    [pPGA,pgaEdges,pgaVect,pgaHat] = pGMP_MC(pgam,mIn,rIn,da,nrnd_pgx);
                    [pPGV,pgvEdges,pgvVect,pgvHat] = pGMP_MC(pgvm,mIn,rIn,da,nrnd_pgx);
                    iMMI                           = pgx2immi(10.^pgaVect,10.^pgvVect);
                    pimmi                          = sampleVect2pdf(iMMI.immi,dimmi);
                    %clf; plot(pgvEdges(1:end-1)+da/2,pPGV); hold on; plot(pgvHat.val,pgvHat.p,'or')
                    
                    % Save PGx distributions
                    ssrte.pythia_ss.pga.prct{it} = prctile(pgaVect,[1,5,16,50,84,95,99]);
                    ssrte.pythia_ss.pga.hat (it) = pgaHat.val;
                    ssrte.pythia_ss.pga.mean(it) = mean(pgaVect);
                    ssrte.pythia_ss.pga.var (it) = var (pgaVect);
                    ssrte.pythia_ss.pgv.prct{it} = prctile(pgvVect,[1,5,16,50,84,95,99]);
                    ssrte.pythia_ss.pgv.hat (it) = pgvHat.val;
                    ssrte.pythia_ss.pgv.mean(it) = mean(pgvVect);
                    ssrte.pythia_ss.pgv.var (it) = var (pgvVect);
                    
                    % Turn PGx distributions into iMMI distributions
                    ssrte.pythia_ss.mmi.hat (it)  = pimmi.mpe;
                    ssrte.pythia_ss.mmi.mean(it)  = mean(iMMI.immi);
                    ssrte.pythia_ss.mmi.var (it)  = var (iMMI.immi);
                    ssrte.pythia_ss.mmi.prct{it}  = prctile(iMMI.immi,[1,5,16,50,84,95,99]);
                    ssrte.pythia_ss.mmi.pEx(:,it) = get_propOfExceedence_ecdf(iMMI.immi,mmiVectPex);
                end
                
                
                if ismember('pythian_ss',rteMethodList)
                    if o.verbose; fprintf(1,'\n pGMP_MC for Pythia2...'); end
                    %nrnd_pgx=1e5; da=2e-2;
                    [pPGA,pgaEdges,pgaVect,pgaHat] = pGMP_MC_ngaw2('PGA',mIn,rIn,da,nrnd_pgx,vs30target,ngaParams);
                    [pPGV,pgvEdges,pgvVect,pgvHat] = pGMP_MC_ngaw2('PGV',mIn,rIn,da,nrnd_pgx,vs30target,ngaParams);
                    iMMI                           = pgx2immi(10.^pgaVect,10.^pgvVect);
                    pimmi                          = sampleVect2pdf(iMMI.immi,dimmi);
                    %pgaCenters      = pgaEdges(1:end-1)+da/2; clf; plot(pgaCenters,pPGA)
                    %plot(pgvEdges(1:end-1)+da/2,pPGV,'-k'); hold on; plot(pgvHat.val,pgvHat.p,'dr')
                    
                    % Save PGx distributions
                    ssrte.pythian_ss.pga.prct{it} = prctile(pgaVect,[1,5,16,50,84,95,99]);
                    ssrte.pythian_ss.pga.hat (it) = pgaHat.val;
                    ssrte.pythian_ss.pga.mean(it) = mean(pgaVect);
                    ssrte.pythian_ss.pga.var (it) = var(pgaVect);
                    ssrte.pythian_ss.pgv.prct{it} = prctile(pgvVect,[1,5,16,50,84,95,99]);
                    ssrte.pythian_ss.pgv.hat (it) = pgvHat.val;
                    ssrte.pythian_ss.pgv.mean(it) = mean(pgvVect);
                    ssrte.pythian_ss.pgv.var (it) = var(pgvVect);
                    
                    % Turn PGx distributions into iMMI distributions
                    ssrte.pythian_ss.mmi.hat (it)  = pimmi.mpe;
                    ssrte.pythian_ss.mmi.mean(it)  = mean(iMMI.immi);
                    ssrte.pythian_ss.mmi.var (it)  = var (iMMI.immi);
                    ssrte.pythian_ss.mmi.prct{it}  = prctile(iMMI.immi,[1,5,16,50,84,95,99]);
                    ssrte.pythian_ss.mmi.pEx(:,it) = get_propOfExceedence_ecdf(iMMI.immi,mmiVectPex);
                end
                
                
                % Extract Pd-observation for multi-station magnitude comps; use a maximum of 4sec to
                % measure Pd, and truncate time window if S-wave arrives earlier
                nsnp_4sec          = round(4/snpLength);
                lastSnippetBeforeS = floor((zTarget.px.s.t-zTarget.px.p.t)/snpLength);
                nsnp_max_pd        = min([nsnp_4sec, lastSnippetBeforeS]);
                
                if nsnpavail<=nsnp_max_pd; nsnp_pd = nsnpavail;
                else                       nsnp_pd = nsnp_max_pd;
                end
                
                % TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP
                if nsnp_pd<0; 
                    nsnp_pd=4;
                    fprintf(1,'this should no longer happen after NGA-bug has been corrected.\n')
                end
                %pdi                                  = zTarget.var8{1}(nsnp_pd); % --> max. pd from all 3 comps
                pdi                                  = zTarget.pgd.tszen{1}(nsnp_pd); % --> max. pd from all 3 comps
                rte.ka13.pd.val(iiCurrentTrigger,it) = pdi;
                
                % Make single station ka13-type GMP
                if ismember('ka13n_ss',rteMethodList)
                    fprintf('rIn.var is computed using sample variance of log-normal distribution. Make sure this is correctly translated into location uncertainty for GMPE... not yet done. do it. now.\n')
                    pause
                    
                    rHatPd  = 10^gbaMap.r.mpe;
                    mHatPd  = 1.23*log10(pdi) +1.38*log10(rHatPd) +5.39;
                    mIn.mu  = mHatPd;
                    mIn.var = .4^2;
                    rIn.mu  = rHatPd;
                    rIn.var = gbaMap.r.hat2.std^2;
                    
                    [pPGA,pgaEdges,pgaVect,pgaHat] = pGMP_MC_ngaw2('PGA',mIn,rIn,da,nrnd_pgx,vs30target,ngaParams);
                    [pPGV,pgvEdges,pgvVect,pgvHat] = pGMP_MC_ngaw2('PGV',mIn,rIn,da,nrnd_pgx,vs30target,ngaParams);
                    iMMI                           = pgx2immi(10.^pgaVect,10.^pgvVect);
                    pimmi                          = sampleVect2pdf(iMMI.immi,dimmi);
                    
                    ssrte.ka13n_ss.m.hat   (it)  = mHatPd;
                    ssrte.ka13n_ss.mmi.hat (it)  = pimmi.mpe;
                    ssrte.ka13n_ss.mmi.mean(it)  = mean(iMMI.immi);
                    ssrte.ka13n_ss.mmi.var (it)  = var (iMMI.immi);
                    ssrte.ka13n_ss.mmi.prct{it}  = prctile(iMMI.immi,[1,5,16,50,84,95,99]);
                    ssrte.ka13n_ss.mmi.pEx(:,it) = get_propOfExceedence_ecdf(iMMI.immi,mmiVectPex);
                end
                
                
                % Save ssrte to array
                ssrteArray{iiCurrentTrigger} = ssrte;
                clear ssrte

                %                 if o.writeFinderFiles
                %                     paObs                             = zTarget.var4{1}(nsnpavail); % var4 contains max amp from 3 components
                %                     finderObsMat (iiCurrentTrigger,3) = log10(paObs*100);   % Finder needs log10(pa) in cm/s/s
                %                     logPaPred                         = gbaMap.pga.mean;    % pga.mean already comes as log10
                %                     finderPredMat(iiCurrentTrigger,3) = logPaPred+2;        % add 2 to translate from m/s/s to cm/s/s
                %                 end 
            end
            
            
            %% Offsite/multi-station estimates for all stations
            if ntrig~=0
                
                
                if o.verbose; fprintf(1,'\n Multi-station location estimates...'); end
                
                % LOCATION ESTIMATES ..................................
                % Perturb Epicenter Location
                if     ntrig==1; stdLoc = 30;
                elseif ntrig==2; stdLoc = 20;
                elseif ntrig==3; stdLoc = 10;
                elseif ntrig==4; stdLoc = 8;
                elseif ntrig==5; stdLoc = 6;
                elseif ntrig>=6; stdLoc = 5;
                end
                
                dkm      = stdLoc/sqrt(2);
                km2lon   = 1/lon2km(1,eqLat);
                sigmaLat = km2deg(dkm);
                sigmaLon = dkm*km2lon;
                dlat     = sigmaLat*randn(1,1);
                dlon     = sigmaLon*randn(1,1);
                eqLatHat = eqLat+dlat;
                eqLonHat = eqLon+dlon;
                
                rte.loc.rand.std(it)     = stdLoc;
                rte.loc.rand.lat.hat(it) = eqLatHat;
                rte.loc.rand.lon.hat(it) = eqLonHat;
                
                rHypHat = zeros(ntargets,1);
                rEpiHat = zeros(ntargets,1);
                for itarget = 1:ntargets
                    
                    stLat = zTargetList.station.lat(itarget);
                    stLon = zTargetList.station.lon(itarget);
                    stAlt = zTargetList.station.alt(itarget);
                    
                    rHypHat(itarget) = hypoDistance(eqLatHat,eqLonHat,zDefault,stLat,stLon,stAlt)/1000;
                    rEpiHat(itarget) = hypoDistance(eqLatHat,eqLonHat,0       ,stLat,stLon,0    )/1000;
                end
                % Compute location based on obs
                % from Lomax et al., 2007
                % locHat1 = get_arrival_time_loc(zTargetList,hasTrig,stationList);
                % plot_location_map(zTargetList,hasTrig,locHat1,stationList,optsMapPlot)
                % locHat2 = get_GbA_loc(inputs);
                % [Nyad] = px2nyad(iSt_wPx,iSt_woPx,px_times,t_run);
                
                
                
                % MAGNITUDE ESTIMATES .................................
                if ismember('gba',rteMethodList)
                    if o.verbose; fprintf(1,'\n Multi-station magnitude estimates...'); end
                    % 1. GbA multi-station magnitude estimates: get magnitude marginal posteriors from
                    % all triggers with sufficient available data; combine and maximise them.
                    % Note: tavail has dimensions of idxTargets, and idxTargets are rel. to zList, in which llh(m) functions are stored
                    
                    % Strategy 1: use llh functions from all triggered stations
                    iiHasData   = (tavail>0);
                    iiAintEmpty = cellfun(@(x) ~isempty(x), lmmi);
                    iiUseMe     = iiHasData &iiAintEmpty;
                    lmmiMat     = cell2mat(lmmi(iiUseMe));
                    lmm         = sum(lmmiMat,1);
                    [~,imax]    = max(lmm);
                    mMpe        = mm(imax);
                    
                    % Find variance of m-pdf
                    lmmc     = lmm./abs(min(lmm(lmm~=-inf)));
                    Lmm      = 10.^lmmc;
                    Lmmn     = Lmm/sum(Lmm);
                    musig0   = [mMpe; 0.5];
                    [musigHat,~,flg] = fminsearch(@(musig) lsq_normDist(musig,mm,Lmmn),musig0);
                    % CHECK: yfit  = normpdf(mm,musigHat(1),musigHat(2)); yfitn = yfit/sum(yfit);
                    %clf; hold on; plot(mm,Lmmn,'-xk'); plot(mm,yfitn,'-dr')
                    
                    rte.gba.m.hat(it) = mMpe;
                    rte.gba.m.std(it) = musigHat(2);
                    
                    
                    
                    % Strategy 2: if eq could be big, use llh functions only from stations with minimum data availability
                    it_last=it-1; if it_last<1; it_last=1; end
                    mHat1_last = rte.gba.m.hat (it_last);
                    mHat2_last = rte.gba.m.hat2(it_last);
                    
                    if mHat1_last<mminGba2 ||mHat2_last<mminGba2; mMpe_all = mMpe;
                    else
                        iiHasData   = (tavail>tavailMin);
                        iiAintEmpty = cellfun(@(x) ~isempty(x), lmmi);
                        iiUseMe     = iiHasData &iiAintEmpty;
                        
                        % If no stations have tavailMin, use the one station with the most data
                        if sum(iiUseMe)==0; [val,iiMax] = max(tavail);
                            iiUseMe     = iiMax;  % check: tavail(iiMax)==max(tavail)
                        end
                        lmmiMat  = cell2mat(lmmi(iiUseMe));
                        lmm      = sum(lmmiMat,1);
                        [~,imax] = max(lmm);
                        mMpe_all = mm(imax);
                    end
                    
                    % Find variance of m-pdf
                    lmmc     = lmm./abs(min(lmm(lmm~=-inf)));
                    Lmm      = 10.^lmmc;
                    Lmmn     = Lmm/sum(Lmm);
                    musig0   = [mMpe; 0.5];
                    [musigHat,~,flg] = fminsearch(@(musig) lsq_normDist(musig,mm,Lmmn),musig0);
                    if flg==0;
                        1+1;
                    end
                    rte.gba.m.hat2(it) = mMpe_all;
                    rte.gba.m.std2(it) = musigHat(2);
                end
                
                
                % 2. Pd-magnitude estimates from Kuyuk & Allen, 2013, GRL ("ka13")
                %    use 4s since P-trigger and avoid S-phases
                pdVal        = rte.ka13.pd.val(:,it)*100; % [cm]
                hasTrig      = pdVal~=0;
                pdVal        = pdVal  (hasTrig);
                rEpiHat_trig = rEpiHat(hasTrig);
                mPdVect      = 1.23*log10(pdVal) +1.38*log10(rEpiHat_trig) +5.39;
                mPdHat       = mean(mPdVect);
                rte.ka13.m.hat (it) = mPdHat;
                rte.ka13.m.val{it}  = mPdVect;
                rte.ka13.m.std(it)  = 0.31;
                
                
                % GROUND MOTION PREDICTIONS ...........................
                % GMP at all sites, inlcuding non-triggered ones, based
                % on the multi-station estimates for m & r
                
                fprintf(1,' Off-site GMPs...');
                
                %profile on 
                for itarget = 1:ntargets
                    
                    ssrte = ssrteArray{itarget};
                    %ssrte = zList.var.v1{idxTargets(itarget)};
                    vs30target = zTargetList.station.vs30(itarget);
                    
                    % Distance estimate and distance error for this site
                    rEpiHat1 = rEpiHat(itarget);
                    rHypHat1 = rHypHat(itarget);
                    dr_hyp   = rHypHat1-zTargetList.dist.hyp(itarget);
                    dr_epi   = rEpiHat1-zTargetList.dist.epi(itarget);
                    rte.loc.rand.dr_hyp(itarget,it) = dr_hyp;
                    rte.loc.rand.dr_epi(itarget,it) = dr_epi;
                    
                    
                    % GMP for this site ...
                    % ... using GbA
                    mIn.mu  = mMpe;
                    mIn.var = rte.gba.m.std(it)^2;
                    rIn.mu  = rHypHat1;
                    rIn.var = stdLoc^2;     % IS IT OK TO SIMPLY RAISE THE LOG-VARIANCE?

                    if ismember('gba',rteMethodList)
                        [pPGA,pgaEdges,pgaVect,pgaHat]  = pGMP_MC(pgam,mIn,rIn,da,nrnd_pgx);
                        [pPGV,pgvEdges,pgvVect,pgvHat]  = pGMP_MC(pgvm,mIn,rIn,da,nrnd_pgx);
                        iMMI                            = pgx2immi(10.^pgaVect,10.^pgvVect);
                        pimmi                           = sampleVect2pdf(iMMI.immi,dimmi);
                        
                        ssrte.gba.pga.hat  (it)  = pgaHat.val;
                        ssrte.gba.pgv.hat  (it)  = pgvHat.val;
                        ssrte.gba.m.hat    (it)  = mMpe;
                        ssrte.gba.m.hat2   (it)  = mMpe_all;
                        ssrte.gba.mmi.hat (it)  = pimmi.mpe;
                        ssrte.gba.mmi.mean(it)  = mean(iMMI.immi);
                        ssrte.gba.mmi.var (it)  = var (iMMI.immi);
                        ssrte.gba.mmi.prct{it}  = prctile(iMMI.immi,[1,5,16,50,84,95,99]);
                        ssrte.gba.mmi.pEx(:,it) = get_propOfExceedence_ecdf(iMMI.immi,mmiVectPex);
                        
                        % Save multi-station GbA magnitudes estimates to single-site rte
                        ssrte.gba.m.hat2(it) = mMpe_all;
                        ssrte.gba.m.hat (it) = mMpe;
                    end
                    
                    % Also store point-estimate immiHat for prediction with NGA-models
                    if ismember('gban',rteMethodList)
                        fprintf('Should be using GbA estimate for rIn here... not yet done. do it. now.\n')
                        pause
                        
                        [pPGA,pgaEdges,pgaVect,pgaHat] = pGMP_MC_ngaw2('PGA',mIn,rIn,da,nrnd_pgx,vs30target,ngaParams);
                        [pPGV,pgvEdges,pgvVect,pgvHat] = pGMP_MC_ngaw2('PGV',mIn,rIn,da,nrnd_pgx,vs30target,ngaParams);
                        iMMI                           = pgx2immi(10.^pgaVect,10.^pgvVect);
                        pimmi                          = sampleVect2pdf(iMMI.immi,dimmi);
                        
                        ssrte.gban.pga.hat   (it) = pgaHat.val;
                        ssrte.gban.pgv.hat   (it) = pgvHat.val;
                        ssrte.gban.m.hat     (it)  = mMpe;
                        ssrte.gban.m.hat2    (it)  = mMpe_all;
                        ssrte.gban.mmi.hat  (it) = pimmi.mpe;
                        ssrte.gban.mmi.mean (it) = mean(iMMI.immi);
                        ssrte.gban.mmi.var  (it) = var (iMMI.immi);
                        ssrte.gban.mmi.prct {it} = prctile(iMMI.immi,[1,5,16,50,84,95,99]);
                        ssrte.gban.mmi.pEx(:,it) = get_propOfExceedence_ecdf(iMMI.immi,mmiVectPex);
                    end
                    
                    % Also store point-estimate immiHat for prediction with NGA-models
                    if ismember('gbanv',rteMethodList)
                        fprintf('Should be using GbA estimate for rIn here... not yet done. do it. now.\n')
                        pause
                        
                        [pPGA,pgaEdges,pgaVect,pgaHat] = pGMP_MC_ngaw2('PGA',mIn,rIn,da,nrnd_pgx,vs30target,ngaParams);
                        [pPGV,pgvEdges,pgvVect,pgvHat] = pGMP_MC_ngaw2('PGV',mIn,rIn,da,nrnd_pgx,vs30target,ngaParams);
                        iMMI                           = pgx2immi(10.^pgaVect,10.^pgvVect);
                        pimmi                          = sampleVect2pdf(iMMI.vmmi,dimmi);
                        
                        ssrte.gban.pga.hat  (it) = pgaHat.val;
                        ssrte.gban.pgv.hat  (it) = pgvHat.val;
                        ssrte.gban.m.hat    (it) = mMpe;
                        ssrte.gban.m.hat2   (it) = mMpe_all;
                        ssrte.gban.mmi.hat  (it) = pimmi.mpe;
                        ssrte.gban.mmi.mean (it) = mean(iMMI.vmmi);
                        ssrte.gban.mmi.var  (it) = var (iMMI.vmmi);
                        ssrte.gban.mmi.prct {it} = prctile(iMMI.vmmi,[1,5,16,50,84,95,99]);
                        ssrte.gban.mmi.pEx(:,it) = get_propOfExceedence_ecdf(iMMI.vmmi,mmiVectPex);
                    end
                    
                    
                    
                    % ... using m(Pd)
                    mIn.mu  = mPdHat;
                    mIn.var = rte.ka13.m.std(it)^2;
                    rIn.mu  = rHypHat1;
                    rIn.var = stdLoc^2;
                    
                    if ismember('ka13',rteMethodList)
                        [pPGA,pgaEdges,pgaVect,pgaHat] = pGMP_MC(pgam,mIn,rIn,da,nrnd_pgx);
                        [pPGV,pgvEdges,pgvVect,pgvHat] = pGMP_MC(pgvm,mIn,rIn,da,nrnd_pgx);
                        iMMI                           = pgx2immi(10.^pgaVect,10.^pgvVect);
                        pimmi                          = sampleVect2pdf(iMMI.immi,dimmi);
                        
                        ssrte.ka13.pga.hat  (it) = pgaHat.val;
                        ssrte.ka13.pgv.hat  (it) = pgvHat.val;
                        ssrte.ka13.m.hat    (it) = mPdHat;
                        ssrte.ka13.mmi.hat (it)  = pimmi.mpe;
                        ssrte.ka13.mmi.mean(it)  = mean(iMMI.immi);
                        ssrte.ka13.mmi.var (it)  = var (iMMI.immi);
                        ssrte.ka13.mmi.pEx(:,it) = get_propOfExceedence_ecdf(iMMI.immi,mmiVectPex);
                        ssrte.ka13.mmi.prct{it}  = prctile(iMMI.immi,[1,5,16,50,84,95,99]);
                        ca = pgam.coeffs; cv = pgvm.coeffs;
                        pgaHat                    = pgam_rock.fgmp(rHypHat1,mPdHat,ca(1),ca(2),ca(3),ca(4),ca(5),ca(6));
                        pgvHat                    = pgvm_rock.fgmp(rHypHat1,mPdHat,cv(1),cv(2),cv(3),cv(4),cv(5),cv(6));
                    end
                    
                    if ismember('ka13n',rteMethodList)
                        [pPGA,pgaEdges,pgaVect,pgaHat] = pGMP_MC_ngaw2('PGA',mIn,rIn,da,nrnd_pgx,vs30target,ngaParams);
                        [pPGV,pgvEdges,pgvVect,pgvHat] = pGMP_MC_ngaw2('PGV',mIn,rIn,da,nrnd_pgx,vs30target,ngaParams);
                        iMMI                           = pgx2immi(10.^pgaVect,10.^pgvVect);
                        pimmi                          = sampleVect2pdf(iMMI.immi,dimmi);
                        
                        ssrte.ka13n.pga.hat (it)  = pgaHat.val;
                        ssrte.ka13n.pgv.hat (it)  = pgvHat.val;
                        ssrte.ka13n.m.hat   (it)  = mPdHat;
                        ssrte.ka13n.mmi.hat (it)  = pimmi.mpe;
                        ssrte.ka13n.mmi.mean(it)  = mean(iMMI.immi);
                        ssrte.ka13n.mmi.var (it)  = var (iMMI.immi);
                        ssrte.ka13n.mmi.pEx(:,it) = get_propOfExceedence_ecdf(iMMI.immi,mmiVectPex);
                        ssrte.ka13n.mmi.prct{it}  = prctile(iMMI.immi,[1,5,16,50,84,95,99]);
                        %ssrte.ka13.mmi.hat3(it)        = prediction wo/ error propagation;
                    end
                    
                    if ismember('ka13nv',rteMethodList)
                        [pPGA,pgaEdges,pgaVect,pgaHat] = pGMP_MC_cuh07m2(cuh.pga.rhyp,mIn,rIn,da,nrnd_pgx);
                        [pPGV,pgvEdges,pgvVect,pgvHat] = pGMP_MC_cuh07m2(cuh.pgv.rhyp,mIn,rIn,da,nrnd_pgx);
                        %[pPGA,pgaEdges,pgaVect,pgaHat] = pGMP_MC_ngaw2('PGA',mIn,rIn,da,nrnd_pgx,vs30target,ngaParams);
                        %[pPGV,pgvEdges,pgvVect,pgvHat] = pGMP_MC_ngaw2('PGV',mIn,rIn,da,nrnd_pgx,vs30target,ngaParams);
                        iMMI                           = pgx2immi(10.^pgaVect,10.^pgvVect);
                        pimmi                          = sampleVect2pdf(iMMI.vmmi,dimmi);
                        
                        ssrte.ka13nv.pga.hat (it)  = pgaHat.val;
                        ssrte.ka13nv.pgv.hat (it)  = pgvHat.val;
                        ssrte.ka13nv.m.hat   (it)  = mPdHat;
                        ssrte.ka13nv.mmi.hat (it)  = pimmi.mpe;
                        ssrte.ka13nv.mmi.mean(it)  = mean(iMMI.vmmi);
                        ssrte.ka13nv.mmi.var (it)  = var (iMMI.vmmi);
                        ssrte.ka13nv.mmi.pEx(:,it) = get_propOfExceedence_ecdf(iMMI.vmmi,mmiVectPex);
                        %ssrte.ka13nv.mmi.prct{it}  = prctile(iMMI.vmmi,[1,5,16,50,84,95,99]);
                    end
                    
                    if ismember('ka13nva',rteMethodList)
                    
                        pimmi                       = sampleVect2pdf(iMMI.immi,dimmi);
                        ssrte.ka13nva.pga.hat (it)  = pgaHat.val;
                        ssrte.ka13nva.pgv.hat (it)  = pgvHat.val;
                        ssrte.ka13nva.m.hat   (it)  = mPdHat;
                        ssrte.ka13nva.mmi.hat (it)  = pimmi.mpe;
                        ssrte.ka13nva.mmi.mean(it)  = mean(iMMI.immi);
                        ssrte.ka13nva.mmi.var (it)  = var (iMMI.immi);
                        ssrte.ka13nva.mmi.pEx(:,it) = get_propOfExceedence_ecdf(iMMI.immi,mmiVectPex);
                        %ssrte.ka13nva.mmi.prct{it}  = prctile(iMMI.immi,[1,5,16,50,84,95,99]);
                    end
                    
                    % ... using catalog values: pythian
                    mIn.mu  = m;
                    mIn.var = 0;
                    rIn.mu  = zTargetList.dist.flt(itarget);
                    rIn.var = 0;

                    if ismember('pythia',rteMethodList)
                        [pPGA,pgaEdges,pgaVect,pgaHat] = pGMP_MC(pgam,mIn,rIn,da,nrnd_pgx);
                        [pPGV,pgvEdges,pgvVect,pgvHat] = pGMP_MC(pgvm,mIn,rIn,da,nrnd_pgx);
                        iMMI                           = pgx2immi(10.^pgaVect,10.^pgvVect);
                        pimmi                          = sampleVect2pdf(iMMI.immi,dimmi);
                        
                        ssrte.pythia.pga.hat (it)  = pgaHat.val;
                        ssrte.pythia.pgv.hat (it)  = pgvHat.val;
                        ssrte.pythia.mmi.hat (it)  = pimmi.mpe;
                        ssrte.pythia.mmi.mean(it)  = mean(iMMI.immi);
                        ssrte.pythia.mmi.var (it)  = var (iMMI.immi);
                        ssrte.pythia.mmi.pEx(:,it) = get_propOfExceedence_ecdf(iMMI.immi,mmiVectPex);
                        ssrte.pythia.mmi.prct{it}  = prctile(iMMI.immi,[1,5,16,50,84,95,99]);
                        % ca = pgam.coeffs; cv = pgvm.coeffs;
                        %pgaHat                    = pgam_rock.fgmp(rHypHat1,mPdHat,ca(1),ca(2),ca(3),ca(4),ca(5),ca(6));
                        %pgvHat                    = pgvm_rock.fgmp(rHypHat1,mPdHat,cv(1),cv(2),cv(3),cv(4),cv(5),cv(6));
                    end
                    
                    % Also store point-estimate immiHat for prediction with NGA-models
                    if ismember('pythian',rteMethodList)
                        [pPGA,pgaEdges,pgaVect,pgaHat] = pGMP_MC_ngaw2('PGA',mIn,rIn,da,nrnd_pgx,vs30target,ngaParams);
                        [pPGV,pgvEdges,pgvVect,pgvHat] = pGMP_MC_ngaw2('PGV',mIn,rIn,da,nrnd_pgx,vs30target,ngaParams);
                        iMMI                           = pgx2immi(10.^pgaVect,10.^pgvVect);
                        pimmi                          = sampleVect2pdf(iMMI.immi,dimmi);
                        
                        ssrte.pythian.pga.hat (it)  = pgaHat.val;
                        ssrte.pythian.pgv.hat (it)  = pgvHat.val;
                        ssrte.pythian.mmi.hat (it)  = pimmi.mpe;
                        ssrte.pythian.mmi.mean(it)  = mean(iMMI.immi);
                        ssrte.pythian.mmi.var (it)  = var (iMMI.immi);
                        ssrte.pythian.mmi.pEx(:,it) = get_propOfExceedence_ecdf(iMMI.immi,mmiVectPex);
                        ssrte.pythian.mmi.prct{it}  = prctile(iMMI.immi,[1,5,16,50,84,95,99]);
                    end
                    
                    if ismember('pythianv',rteMethodList)
                        [pPGA,pgaEdges,pgaVect,pgaHat] = pGMP_MC_cuh07m2(cuh.pga.rflt,mIn,rIn,da,nrnd_pgx);
                        [pPGV,pgvEdges,pgvVect,pgvHat] = pGMP_MC_cuh07m2(cuh.pgv.rflt,mIn,rIn,da,nrnd_pgx);
                        %[pPGA,pgaEdges,pgaVect,pgaHat] = pGMP_MC_ngaw2('PGA',mIn,rIn,da,nrnd_pgx,vs30target,ngaParams);
                        %[pPGV,pgvEdges,pgvVect,pgvHat] = pGMP_MC_ngaw2('PGV',mIn,rIn,da,nrnd_pgx,vs30target,ngaParams);
                        iMMI                           = pgx2immi(10.^pgaVect,10.^pgvVect);
                        pimmi                          = sampleVect2pdf(iMMI.vmmi,dimmi);
                        
                        ssrte.pythianv.pga.hat (it)  = pgaHat.val;
                        ssrte.pythianv.pgv.hat (it)  = pgvHat.val;
                        ssrte.pythianv.mmi.hat (it)  = pimmi.mpe;
                        ssrte.pythianv.mmi.mean(it)  = mean(iMMI.vmmi);
                        ssrte.pythianv.mmi.var (it)  = var (iMMI.vmmi);
                        ssrte.pythianv.mmi.pEx(:,it) = get_propOfExceedence_ecdf(iMMI.vmmi,mmiVectPex);
                        %ssrte.pythianv.mmi.prct{it}  = prctile(iMMI.vmmi,[1,5,16,50,84,95,99]);
                    end
                    
                    if ismember('pythianva',rteMethodList)
                        pimmi                         = sampleVect2pdf(iMMI.immi,dimmi);
                        ssrte.pythianva.pga.hat (it)  = pgaHat.val;
                        ssrte.pythianva.pgv.hat (it)  = pgvHat.val;
                        ssrte.pythianva.mmi.hat (it)  = pimmi.mpe;
                        ssrte.pythianva.mmi.mean(it)  = mean(iMMI.immi);
                        ssrte.pythianva.mmi.var (it)  = var (iMMI.immi);
                        ssrte.pythianva.mmi.pEx(:,it) = get_propOfExceedence_ecdf(iMMI.immi,mmiVectPex);
                        %ssrte.pythianva.mmi.prct{it}  = prctile(iMMI.immi,[1,5,16,50,84,95,99]);
                    end
                    
                    ssrteArray{itarget}=ssrte ;
                    %zList.var.v1{idxTargets(itarget)} = ssrte;
                end
                %profile viewer
                fprintf(1,' done.');
                
            else
                mMpe     = 0;
                mMpe_all = 0;
            end
            
            
            % Write finder file
            %             if o.writeFinderFiles
            %                 finderObsSnpFullName  = sprintf('%sdata_%i',finderObsFullName ,it-1);
            %                 finderPredSnpFullName = sprintf('%sdata_%i',finderPredFullName,it-1);
            %                 dlmwrite(finderObsSnpFullName ,finderObsMat ,'delimiter',' ','precision',9)
            %                 dlmwrite(finderPredSnpFullName,finderPredMat,'delimiter',' ','precision',9)
            %             end
            if o.verbose; fprintf(1,'\n done.\n'); end
        end
        
        % Evaluate prediction performance at all sites. For each site extract ssrte, compute 
        % GMPPerformance, write it into ssrte & return updated version of ssrte; clean out 
        % empty fields to save space; overwrite ssrte entry in zList
        fprintf(1,'UNCOMMENT ME!\n')
%         fprintf(1,sprintf('\nEvaluating predictions for all %i records & saving ssrte in zList.var.v1...',ntargets))
%         for itarget = 1:ntargets
%             
%             %gmpOpts.plotGMP=0;
%             print_iteration_numbers(itarget,ntargets,'ones')
%             idxThisCorec             = idxTargets(itarget);
%             ssrte                    = ssrteArray{itarget};
%             [ssrte_updated,~]        = evaluate_gmp_1site(ssrte,zList,idxThisCorec,rteMethodList,immiThreshold,gmpOpts);
%             [ssrte_clean]            = clean_out_ssrte(ssrte_updated); % Clean out ssrte to save space
%             zList.var.v1{idxThisCorec} = ssrte_clean;
%         end
        
        eqs.rte{ieq} = rte;
        clear rte
        
        % Plot a short performance summary for each event; was done for AGU Poster 2014 on Sth Napa event
        % if o.plot_all_traces_of_eq; hf = plot_singleEq_rte_summary_NOTYETRUNNING; end
        
        
        % Write GbA estimates and inputs to ascii files
        if o.writeGbaAsciiOut
            write_gba_out_to_ascii(gbaAscii,gbaAsciiPath,zTargetList)
        end

        
        % Temporary Save
        if ismember(ct_eqProc,tens) &o.saveOut
            tmpOutFullName = sprintf('%sout_tr%i_pid%i_ieq%i',projectDirFullName,nz,pid,ieq);
            fprintf(1,sprintf('Saving temporary zList & eqs structure to \n\t%s ...',tmpOutFullName))
            
            % Save some space by saving only trList- and eqs-entries of large magnitude events
            % 1. TraceLists
            clear out
            out.zList = zList.selectSubList(find(zList.eq.m>mminTarget));
            %out.hList = hList.selectSubList(find(hList.eq.m>mminTarget));
            
            rmFieldList = {'fb.amax'; 'fb.cav'; 'var.v2'; 'var.v3'; 'var.v4'; 'var.v5'; 'var.v6'; 'var.v7'; 'var.v8'; ...
                'pga.p'; 'pgv.p'; 'pgd.p'; 'pga.s'; 'pgv.s'; 'pgd.s'; 'pga.nb'; 'pgv.nb'; 'pgd.nb';'noise.nb';'px.allPx'};
            overwrite_unneeded_fields(out.zList,rmFieldList);

            % 2. eqs-structure
            eqs_out = compile_eqs_structure(out.zList,[]);   % For 'traceId' & 'eqIdx'
            eqs_out = merge_double_eqs_entries(eqs_out,'3031111','LANDERS_1992');
            eqs_out = merge_double_eqs_entries(eqs_out,'9108652','HECTOR_1999' );
            eqs_out = merge_double_eqs_entries(eqs_out,'3144585','NORTHR_1994' );
            neqtmp  = numel(eqs_out.m);
            eqs_out.rte = cell(neqtmp,1); 
            for ieq=1:neqtmp
                hasSameName = strcmp(eqs_out.name{ieq},eqs.name);
                hasSameMagn = eqs_out.m(ieq)==eqs.m;
                hasSameLat  = eqs_out.lat(ieq)==eqs.lat;
                hasSameLon  = eqs_out.lon(ieq)==eqs.lon;
                idxMatch    = find(hasSameName &hasSameMagn &hasSameLat &hasSameLon);
                if numel(idxMatch)~=1; fprintf(1,'8UNG: ambiguous eqs-match!\n'); ct_eqs = ct_eqs+1; end
                eqs_out.rte{ieq} = eqs.rte{idxMatch};
            end
            out.zList.prop.eqs         = eqs_out;
            out.zList.prop.stationList = stationList;
            %eqs_out = extract_subeqs(eqs,mminTarget);        % For all fieds except 'traceId' & 'eqIdx'
            %             if ~isequal(eqs_out.m,eqs_tmp.m); error('Something wrong with eqs-structure. Check.'); end
            %             eqs_out.traceId = eqs_tmp.traceId;
            %             eqs_out.eq.idx   = eqs_tmp.eq.idx;
                        % CHECK: ieq = 32; idxList = eqs_out.traceId{ieq}; tmpList.fullName(idxList)
            
            % Save more information
            out.pid         = pid;
            out.snpList     = snpList;
            out.snpLength   = snpLength;
            out.mmiThresholds = mmiThresholdVect;
            out.tms         = tms;              % time vector for multi station inference
            out.rmaxTarget  = rmaxTarget;
            out.mminTarget  = mminTarget;
            out.gmpOpts     = gmpOpts;
            out.zList.prop.mmiVectPex = mmiVectPex;
            tic; save(tmpOutFullName,'out'); toc
            fprintf(1,'done.\n')
        end
    else
    	fprintf(1,sprintf('Eq %s (%i/%i) is not used as target site, m%3.1f\n',eqs.name{ieq},ieq,neq,eqs.m(ieq)))
    end
end % .. for ieq=1:neq


%             .-._   _ _ _ _ _ _ _ _
%  .-''-.__.-'00  '-' ' ' ' ' ' ' ' '-.
% '.___ '    .   .--_'-' '-' '-' _'-' '._
%  V: V 'vv-'   '_   '.       .'  _..' '.'.
%    '=.____.=_.--'   :_.__.__:_   '.   : :
%            (((____.-'        '-.  /   : :
%                              (((-'\ .' /
%                            _____..'  .'
%                           '-._____.-'



%% 5.2.3 Single station inference & Tau_p and Co.
% ===============================================
% if o.singleStation;
%     singleStation_pseudofct;
%     tau_CP_pseudofct.m
% end


%   \____            ____
%    \   \           \   \
%     \SED\_____      \   \
%      \...~-__()______\___\_____________________
%       \         Filter Bank Airlines        ___\
%        \     oo ooooooooooooooooooooo o o  |_O__\_
%         ~~--_________/~~~/________________________)
%              |      /   /()                  |
%              0     /   /()                   0
%                   /___/
%




%% %%%%%%%%%%%%%%%%%
% 6. SAVE RESULTS  %
%%%%%%%%%%%%%%%%%%%%

if o.saveOut
    
    dbstop if error
    
    outFullName = sprintf('%sout_tr%i_pid%i_all%i',projectDirFullName,nz,pid,ieq);
    %outFullName = strcat([projectDirFullName,'out_tr',num2str(nz),'_pid',num2str(pid),'.mat']);
    fprintf(1,['Saving zList & eqs structure to \n\t',outFullName,'\n'])
    
    clear out
    out.zList = zList.selectSubList(find(zList.eq.m>mminTarget));
    %out.hList = hList.selectSubList(find(hList.eq.m>mminTarget));
    
    rmFieldList = {'fb.amax'; 'fb.cav'; 'var.v2'; 'var.v3'; 'var.v4'; 'var.v5'; 'var.v6'; 'var.v7'; 'var.v8'; ...
        'pga.p'; 'pgv.p'; 'pgd.p'; 'pga.s'; 'pgv.s'; 'pgd.s'; 'pga.nb'; 'pgv.nb'; 'pgd.nb';'noise.nb';'px.allPx'};
    overwrite_unneeded_fields(out.zList,rmFieldList);
    
    % 2. eqs-structure
    eqs_out = compile_eqs_structure(out.zList,[]);   % For 'traceId' & 'eqIdx'
    eqs_out = merge_double_eqs_entries(eqs_out,'3031111','LANDERS_1992');
    eqs_out = merge_double_eqs_entries(eqs_out,'9108652','HECTOR_1999' );
    eqs_out = merge_double_eqs_entries(eqs_out,'3144585','NORTHR_1994' );
    neqtmp  = numel(eqs_out.m);
    eqs_out.rte = cell(neqtmp,1);
    for ieq=1:neqtmp
        hasSameName = strcmp(eqs_out.name{ieq},eqs.name);
        hasSameMagn = eqs_out.m(ieq)==eqs.m;
        hasSameLat  = eqs_out.lat(ieq)==eqs.lat;
        hasSameLon  = eqs_out.lon(ieq)==eqs.lon;
        idxMatch    = find(hasSameName &hasSameMagn &hasSameLat &hasSameLon);
        if numel(idxMatch)~=1; fprintf(1,'8UNG: ambiguous eqs-match!\n'); ct_eqs = ct_eqs+1; end
        eqs_out.rte{ieq} = eqs.rte{idxMatch};
    end
    out.zList.prop.eqs         = eqs_out;
	out.zList.prop.stationList = stationList;
            
    %     eqs_out = extract_subeqs(eqs,mminTarget);        % For all fieds except 'traceId' & 'eqIdx'
    %     eqs_tmp = compile_eqs_structure(out.zList,[]);   % For 'traceId' & 'eqIdx'
    %     if ~isequal(eqs_out.m,eqs_tmp.m); error('Something wrong with eqs-structure. Check.'); end
    %     eqs_out.traceId = eqs_tmp.traceId;
    %     eqs_out.eq.idx   = eqs_tmp.eq.idx;
    %     out.zList.prop.eqs = eqs_out;
    % CHECK: ieq = 32; idxList = eqs_out.traceId{ieq}; tmpList.fullName(idxList)
    
    % Save more information
    out.pid            = pid;
    out.snpList        = snpList;
    out.snpLength      = snpLength;
    out.mmiThresholds  = mmiThresholdVect;
    out.tms            = tms;              % time vector for multi station inference
    out.rmaxTarget     = rmaxTarget;
    out.mminTarget     = mminTarget;
    out.mmin           = mmin;
    out.gmpOpts        = gmpOpts;
    out.zList.prop.mmiVectPex = mmiVectPex;
    
    fprintf(1,['Saving zList & eqs structure to \n\t',outFullName,' ...'])
    save(outFullName,'out')
    fprintf(1,'done.\n')
    opts.plotData = 0;
    summary_trList(out.zList,[]      ,opts);
    
    if isMasterScript
        unix(['cp -rp config/',configFileName,' ',projectDirFullName]);
        
        fprintf(1,['Writing logfile\n'])
        currentDir    = pwd;
        [~,starttime] = unix('date');
        scriptName    = mfilename('fullpath'); 
        scriptName    = strcat([scriptName,'.m']);
        logFileName   = strcat([outFullName,'.log']);
        unix(['echo "This traceList was compiled with the following script on ',starttime,'" >  ' logFileName]);
        unix(['echo " ',scriptName,'" >>  ' logFileName]);
        unix(['echo " " >>  ' logFileName]);
        unix(['echo "  executed in ',currentDir,'" >>  ' logFileName]);
        unix(['echo " " >>  ' logFileName]);
        unix(['echo " Also check the config (*.txt) file in same directory as this log file" >>  ' logFileName]);
        unix(['echo " " >>  ' logFileName]);
        unix(['echo " " >>  ' logFileName]);
        unix(['echo " " >>  ' logFileName]);
        unix(['echo "....................................................................." >>  ' logFileName]);
        unix(['cat ', scriptName, ' >> ', logFileName]);
    end
end
fprintf(1,sprintf('\n\nWHEN SAVING TMP-EQS STRUCTURES, ITS GONE WRONG IN %i CASES\n\n',ct_eqs))
























%% %%%%%%%%%%%%%
% X. APPENDIX  %
%%%%%%%%%%%%%%%%
% % Bivariate prior
% if o.usePrior
%     logPrior_shpere = get_bivarPrior_GRsphere(1,20,o.plotPrior,o.printPrior);
%     logPrior_flat   = get_bivarPrior_GRflat(1,20,o.plotPrior,o.printPrior);
% end

% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% WHY ARE THERE SO MANY MORE TRACES WITH INSUFFICIENT SNR SINCE I HAVE
% INTRODUCED GAPWINDOW IN SNR COMPUTATION...?! 
% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% nszList = zList.selectSubList(find(zList.dist.hyp<25));
% nnsz = numel(nszList.eq.m);
% sum(nszList.snr<snrMin);
% 
% sWindow     = 1;              % Window lenght for measuring signal amps after p-pick [sec]
% gWindow     = .5;             % Window length between noise window and p-pick [sec]
% nWindow     = 1;              % Window lenght for measuring noise before p-pick [sec]
% 
% for i = 1:nnsz
%     fprintf(1,sprintf('%i\n',i))
%     [S,meta] = read_any_trace(nszList.fullName{i},nszList,1);
%     vel      = S.vel;
%     
%     [snr1(i),~,~,~]        = get_snr    (vel,meta.ppxIdx,meta.sr,sWindow,gWindow,nWindow);
%     [snr2(i),~,~,~]        = get_snr_bak(vel,meta.ppxIdx,meta.sr,sWindow,nWindow);
% end


%configFileName  = 'gba_test_i36.txt';
%configFileName  = 'i32_singleStationMLE.txt';
%configFileName  = 'tdpa_i35_0p01_fast.txt';
%configFileName  = 'tdpa_i35_0p01.txt';
%configFileName  = 'tdpa_i35_0p01_repick.txt';
%configFileName  = 'mEst_i35.txt';
%configFileName  = 'tdpa_i34_0p01_dh.txt';
%configFileName  = 'tdpa_i34_0p01.txt';
%configFileName  = 'tdpa_i34.txt';
%configFileName  = 'tauC_i34.txt';
%configFileName  = 'tdpa_i33_southNapa.txt';
%configFileName  = 'i34.txt';
%configFileName  = 'i33_plot.txt';
%configFileName  = 'i33.txt';
%configFileName  = 'fullmoon.txt';
%configFileName  = 'i32_withPriors.txt';
%configFileName  = 'tau_p_and_co.txt';
%configFileName  = 'acausal_i32.txt';
%configFileName  = 'southNapaEq.txt';
%configFileName = 'paramInf_i32_r0_2_2_Xp1_mSat4.txt';
%configFileName = 'paramInf_i32.txt';
%configFileName = 'paramInf_i32_debug.txt';
%configFileName = 'compute_eqs_wise_regCoeff_CPBmod2.txt';

%                     hasSsrte = cellfun(@(x) isstruct(x), zList.var.v1(idxUseMe));     % Some of the records with enough data have been blocked,
%                     aintEmpty = cellfun(@(x) ~isempty(x), zList.var.v1(idxUseMe));     % Some of the records with enough data have been blocked,
%                     idxUseMe = idxUseMe(hasSsrte);                                  % e.g. because of the rmaxTarget criterion; leave them out.
%                         
%                         idxUseMe2 = idxTargets(tavail>tavailMin);
%                     end
%                     if isempty(idxUseMe2); [val,ii]  = max(tavail);   % If no stations have tavailMin, use the one station with the most data           
%                                            idxUseMe2 = idxTargets(ii); % check: tavail(ii)==max(tavail)
%                     end
%                     hasSsrte  = cellfun(@(x) isstruct(x), zList.var.v1(idxUseMe2));
%                     idxUseMe2 = idxUseMe2(hasSsrte);
%                     lmmMat2   = cell2mat(cellfun(@(x) x.gba.m.llh{it}', zList.var.v1(idxUseMe2),'uniformOutput',0));





%% %%%%%%%%%%%%%%%%%%%%%%%%
% FINITE SOURCE DISTANCES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finSrcFullName = '~/Documents/eq/data/faults/allFaults_i38.mat';
% load(finSrcFullName);
% iiLargeEqs = find(eqs.m>=6.5);
% nlargeEqs  = numel(iiLargeEqs);
% hd_bak     = zList.dist.hyp;
% 
% for ieq = 1:nlargeEqs
%     
%     % Read incdices of all records for this eq
%     iiThisEq    = iiLargeEqs(ieq); 
%     trIndexList = eqs.traceId{iiThisEq};
%     ieq
%     eqs.name{iiThisEq}
%     nst         = numel(trIndexList);
%     stLat       = zList.station.lat(trIndexList);
%     stLon       = zList.station.lon(trIndexList);
%     stAlt       = zList.station.alt(trIndexList);
%     hd          = zList.dist.hyp(trIndexList);
%     
%     % Load corresponding finite fault model
%     ifault = find(strcmp(eqs.name{iiThisEq},allFaults.eqname));
% 
%     if ~isempty(ifault) &~isempty(allFaults.segments{ifault}.lat)
%        
%         thisFault = allFaults.segments{ifault};
%         npoints   = numel(thisFault.lat);
%         
%         % Recompute all distances
%         fd = zeros(nst,1);
%         for ist=1:nst
%             fdVect = zeros(npoints,1);
%             for ipt=1:npoints
%                 fdist_m     = hypoDistance(thisFault.lat(ipt),thisFault.lon(ipt), ...
%                     thisFault.z(ipt)*1e3,stLat(ist),stLon(ist),-stAlt(ist));   % Input in [m] & [deg], output in [m]
%                 fdVect(ipt) = fdist_m*1e-3;
%             end
%             fd(ist) = min(fdVect);
%             %hdb_tmp = zList.dist.hyp(trIndexList(ist));
%             zList.dist.flt(trIndexList(ist)) = fd(ist);
%             eList.dist.flt(trIndexList(ist)) = fd(ist);
%             nList.dist.flt(trIndexList(ist)) = fd(ist);
%             hList.dist.flt(trIndexList(ist)) = fd(ist);
%         end
%         
%         %plot_finite_fault_segment(thisFault); hold on; plot3(stLon,stLat,zeros(size(stLon)),'vk','markerFaceColor','y','markerSize',15)
%         %clf; plot(fd,hd,'xk')
%     else
%         fprintf(1,'8UNG: No finite fault data for eq %s in %s\n',eqs.name{iiThisEq},finSrcFullName)
%     end
% end
%dr=hd_bak-zList.dist.hyp; clf; plot(zList.eq.m,dr,'xk')





% 5.1 Preparations
% ================

% if o.useWt
%     % Approximate X = gradient ratio from dlog(Y)/dM = X*dlog(Y)/dlog(R) at
%     % snippet = 6 and freq-band = 4. Used for computing weights with kNN.
%     X = 0.5;    % X=0.5 is the one from theoretical considerations
%     % TASK: Copmute X for all freq bands (for cj optimisation)                 <-- XXXXXXXXXXXXXX
%     % evaluate_and_check_weightX(zList,eList,nList,6);
%     %wNames = {'W1','W2','W3'};
%     %k = 50; k=20;
%     wNames = {'W3'};
%     
%     W = get_weightVectors(wNames,zList,k,X);
% else
%     W.w0 = ones(numel(zList.eq.m),1);
% end
% 
% 
% %%  1.5 
% %   ==========================================
% m = zList.eq.m;
% r = zList.dist.hyp;
% 
% nsnp    = numel(snpList);
% nsnpmax = max(snpList);
% 
% if o.MEDIANnbPGV
%     clusterList = get_median_amps_in_cellClusters(zList,hList,o.plotClusters);
% end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  4. REGERESSION COEFFICIENTS  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
% if o.simpleLin; [A] = get_regression_coeffs(zList,eList,nList,'simpleLin',snippet,W.w0); end
% if o.BJF93;     [B] = get_regression_coeffs(zList,eList,nList,'BJF93',snippet); end    
% % Non-linear  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
% if o.CPB81mod1; [Ca1] = get_regression_coeffs(zList,eList,nList,'CPB81mod1',snippet,W.w3); end
% %if o.CPB81mod2; [Ca2] = get_regression_coeffs(zList,eList,nList,'CPB81mod2',snippet,W.w3); end
% if o.CPB81mod3; [Ca3] = get_regression_coeffs(zList,eList,nList,'CPB81mod3',snippet,W.w3); end
% if o.CUH07mod1; [Cu1] = get_regression_coeffs(zList,eList,nList,'CUH07mod1',snippet,W.w3); end
% if o.CUH07mod2; [Cu2] = get_regression_coeffs(zList,eList,nList,'CUH07mod2',snippet,W.w3); end
% 
% % Reproduce Campbell81  . . . . . . . . . . . . . . . . . . . . . . . . . .   
% if o.reprodCPB81; dummySnippet = 6;
%  [CPB81] = get_regression_coeffs(zList,eList,nList,'CPB81pga',dummySnippet); 
% end
% 
% % [Gb] = get_regression_coeffs(zList,eList,nList,'CUH07unbound',snippet); % Unbounded Cua & Heatonend
% % [H] = get_regression_coeffs(zList,eList,nList,'CPB81',snippet);
% % [I] = get_regression_coeffs(zList,eList,nList,'CPB81nosat',snippet);
% % [J] = get_regression_coeffs(zList,eList,nList,'CPB81constrained',snippet);
% 
% % Compare regressin models  . . . . . . . . . . . . . . . . . . . . . . . .   
% if o.cfRegModels
%     hbic    = cfRegressionModels({Ca1;Ca2;Ca3;Ca4;Cu1;Cu2;B},'bic');
%     hsigma2 = cfRegressionModels({Gmod;Hmod},'sigma2');
% end
% 
% %  4.2. Compute eventwise reg-coeffs for cross validation
% if o.compRegCoeffs; compute_eqWise_regCoeffs_pseudo; end

% Load event-specific regression coefficients
% if o.CPB81mod2
%     snptxt            = sprintf('%d',snpList);
%     coeffFileFullName = strcat([outDirFullName,'eqCoeff_CPB81mod2_w3_k',num2str(k), ...
%         '_',num2str(nz),'traces_snp',snptxt,'.mat']);
%     load(coeffFileFullName)
%     if o.checks; plot_coeffDistributions(Ca2,snippet); end
%     % Test wether coeffs are ok:    --> recompute_reg_coeffs_1eq.m <--
%     % Note: a-coeffs are really low because log(a) is plotted. a~1e-5, i.e.
%     % log(a)~(-11)
% end
% 
% % Create amplitude sorted indices lists
% if o.mle_NPEMLE
%     fprintf(1,'\nComputing lists with amp-sorted indices for all snippets and bands...\n ')
%     fprintf(1,'\n\n\t8ung: this probably does not work anymore after the changes introduced with multiple station inference...\n ')
%     pause
%     stdIdx = cell(nbands,nsnpmax);
%     for isnp = 1:nsnpmax
%         fprintf(1,[num2str(isnp),' .. '])
%         if ismember(isnp,snpList)
%             for iband = 1:nbands
%                 a                  = cellfun(@(x) x(iband,isnp), zList.fb.amax);
%                 [astd,idx]         = sort(a);
%                 stdIdx{iband,isnp} = idx;
%             end
%         end
%     end
%     fprintf(1,'\n... done.\n')
%     
%     % Test if we get the right indices
%     if o.checks
%         check_similar_amps_sorting  (zList,stdIdx,iband,snippet)
%         check_similar_amps_selection(zList,stdIdx,iband,snippet,nsimhalf)
%     end
% end
% 
% 
% % 5.1.1 Choose data weights
% if o.useWt
%     Wti  = W.w3replicate;            % integer weights
%     Wtn  = single(zeros(size(Wti))); % dummy values
%     Wt0  = single(zeros(size(Wti))); % 'weights' for unweighted case
% end 

%% Choose and parameterise GM prediction models for PGA & PGV

% % Cua & Heaton, 2007
% fgmp_CUH = @(r,m,a,b,c1,c2,d,e) a*m - b*(r+(atan(m-5)+1.4).*(c1*exp(c2*(m-5)))) ...
%            - d*log10(r+(atan(m-5)+1.4).*(c1*exp(c2*(m-5)))) + e;
% 
% % Copied from get_coeffs_CH2007.m:
% % a=0.733;    b=7.216e-4;  d=1.48;    c1=1.16;   c2=0.96;    e=-0.4202;    varHat=0.3069^2; % PGA, Horizontal, S-phase, Rock
% % a=0.709;    b=2.3878e-3; d=1.4386;  c1=1.722;  c2=0.9560;  e=-2.4525e-2; varHat=0.3261^2; % PGA, H, S, S
% % a=0.861988; b=5.578e-4;  d=1.36760; c1=0.8386; c2=0.98;    e=-2.58053;   varHat=0.2773^2; % PGV, H, S, R
% % a=0.920352; b=1.136e-3;  d=1.4999;  c1=1.8825; c2=0.99246; e=-2.30099;   varHat=0.3114^2; % PGV, H, S, S
% 
%     
% % PGA [m/s/s]: Published coefficients are in [cm/s/s] --> change e
% %a=0.779; b=2.555e-3; c1=1.352; c2=1.478; d=1.105; e=-0.645-log10(100); varHat = 0.308^2;       % --> my own inversion?
% a=0.733;    b=7.216e-4;  d=1.48;    c1=1.16;   c2=0.96;    e=-0.4202-log10(100);    varHat=0.3069^2;       % PGA, Horizontal, S-phase, Rock
% pgam_rock.fgmp   = fgmp_CUH; 
% pgam_rock.coeffs = [a; b; c1; c2; d; e];
% pgam_rock.varHat = varHat;
% 
% a=0.709;    b=2.3878e-3; d=1.4386;  c1=1.722;  c2=0.9560;  e=-2.4525e-2-log10(100); varHat=0.3261^2; % PGA, H, S-phase, Soil
% pgam_soil.fgmp   = fgmp_CUH; 
% pgam_soil.coeffs = [a; b; c1; c2; d; e];
% pgam_soil.varHat = varHat;
% 
% % PGV [m/s]: Published coefficients
% %a=0.836; b=2.324e-3; c1=1.562; c2=2.423; d=1.054; e=-0.338-log10(100); varHat = 0.312^2;
% %a=0.836; b=2.324e-3; c1=1.562; c2=2.423; d=1.054; e=-0.338; varHat = 0.312^2;   % From Thesis
% a=0.861988; b=5.578e-4; c1=0.8386; c2=0.98; d=1.36760; e=-2.58053-log10(100); varHat=0.2773^2; % horizontal, S-phase, PGV, Rock
% pgvm_rock.fgmp   = fgmp_CUH;
% pgvm_rock.coeffs = [a; b; c1; c2; d; e];
% pgvm_rock.varHat = varHat;
% 
% a=0.920352; b=1.136e-3;  d=1.4999;  c1=1.8825; c2=0.99246; e=-2.30099-log10(100); varHat=0.3114^2; % PGV, H, S-phase, Soil
% pgvm_soil.fgmp   = fgmp_CUH;
% pgvm_soil.coeffs = [a; b; c1; c2; d; e];
% pgvm_soil.varHat = varHat;
% 
% % Try model
% %rVect = 5:5:100;
% %mTmp  = 5.0;
% %logpgv=fgmp_CUH(rVect,mTmp,a,b,c1,c2,d,e);
% %clf; plot(rVect,10.^logpgv,'-xk'); set(gca,'xscale','log','yscale','log'); grid on

                    



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MANUALLY COMBINE TRUNCATED TRACELISTS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
% 
% load('/scratch/memeier/fbout/i35/zLists/tmp/ntr13884_n0p0_1p0_1p5_2p0_2p5_3p0_4p0_5p0_6p0_7p0_8p0_12p0/nsList_5000.mat')
% t1List = nsList;
% clear nsList
% load('/scratch/memeier/fbout/i35/zLists/tmp/trunc2_ntr13884_n0p0_1p0_1p5_2p0_2p5_3p0_4p0_5p0_6p0_7p0_8p0_12p0/nsList.eq.mat')
% load('/scratch/memeier/fbout/i35/zLists/tmp/trunc2_ntr13884_n0p0_1p0_1p5_2p0_2p5_3p0_4p0_5p0_6p0_7p0_8p0_12p0/var.mat')
% t2List = nsList;
% clear nsList
% load('/scratch/memeier/fbout/i35/zLists/tmp/trunc3_ntr13884_n0p0_1p0_1p5_2p0_2p5_3p0_4p0_5p0_6p0_7p0_8p0_12p0/nsList_9000.mat')
% t3List = nsList;
% clear nsList;
% 
% idx1 = 1   :5000;
% idx2 = 7001:13884;
% idx3 = 5001:7000;
% 
% ntot = numel(t3List.eq.m);
% nexp = numel(outVar.expVect);
% V3   = cell        (ntot,nexp);      % Use for PD matrix
% V4   = single(zeros(ntot,nexp));     % Use for DS matrix
% nsList = t1List;
% 
% V3(idx1,:)  = t1List.var3(idx1,:);
% V3(idx2,:)  = t2List.var3(idx2,:);
% V3(idx3,:)  = t3List.var3(idx3,:);
% V4(idx1,:)  = t1List.var4(idx1,:);
% V4(idx2,:)  = t2List.var4(idx2,:);
% V4(idx3,:)  = t3List.var4(idx3,:);
% nsList.var3 = V3;
% nsList.var4 = V4;
% 
% 
% % Compose file names for saving modified nsList and additional
% % information
% expVect        = outVar.expVect;
% outDirName     = sprintf('/scratch/memeier/fbout/i%i/zLists/',iN);
% exponentString = strrep(sprintf('%3.1f_',expVect),'.','p');
% outName        = sprintf('ntr%i_n%s',ntot,exponentString(1:end-1));
% outDirFullName = sprintf('%s%s/',outDirName,outName);
% outFullName    = sprintf('%s%s',outDirFullName,'nsList.eq.mat');
% 
% % Save final version of nsList
% save(outFullName,'nsList')
% outVarFullName = sprintf('%s%s',outDirFullName,'var.mat');
% save(outVarFullName ,'outVar')




% % MANUALLY COMBINE TWO LISTS
% %a: ntr13889_n1p5_2p5_6p0_7p0_8p0_9p0
% %b: ntr13889_n0p0_1p0_2p0_3p0_4p0_5p0
% 
% load('/scratch/memeier/fbout/i35/zLists/ntr13889_n1p5_2p5_6p0_7p0_8p0_9p0/nsList.eq.mat')
% load('/scratch/memeier/fbout/i35/zLists/ntr13889_n1p5_2p5_6p0_7p0_8p0_9p0/var.mat')
% 
% %a_expVect = expVect;
% aList     = nsList;
% a_outVar  = outVar;
% 
% clear nsList outVar
% load('/scratch/memeier/fbout/i35/zLists/ntr13889_n0p0_1p0_2p0_3p0_4p0_5p0/nsList.eq.mat')
% load('/scratch/memeier/fbout/i35/zLists/ntr13889_n0p0_1p0_2p0_3p0_4p0_5p0/var.mat')
%     
% bList    = nsList;
% b_outVar = outVar;
% 
% expVect = [0 1 1.5 2 2.5 3 4 5 6 7 8 9];    % --> combined exponent vector
% aidx    = [3 5 9 10 11 12];                 % target indices for different columns
% bidx    = [1 2 4 6 7 8];
% 
% ntot = numel(expVect);
% ntmp = numel(nsList.eq.m);
% V3   = cell        (ntmp,ntot);      % Use for PD matrix
% V4   = single(zeros(ntmp,ntot));     % Use for DS matrix
% 
% V3(:,aidx)  = aList.var3;
% V3(:,bidx)  = bList.var3;
% V4(:,aidx)  = aList.var4;
% V4(:,bidx)  = bList.var4;
% nsList.var3 = V3;
% nsList.var4 = V4;
% outVar.expVect       = expVect;
% outVar.saveTimeStamp = clock;
% outVar.nrnd          = 100;
% 
% outDirName     = sprintf('/scratch/memeier/fbout/i%i/zLists/',iN);
% exponentString = strrep(sprintf('%3.1f_',expVect),'.','p');
% outName        = sprintf('ntr%i_n%s',ntmp,exponentString(1:end-1));
% outDirFullName = sprintf('%s%s/',outDirName,outName);
% outFullName    = sprintf('%s%s',outDirFullName,'nsList.eq.mat');
% outVarFullName = sprintf('%s%s',outDirFullName,'var.mat');
% save(outFullName,'nsList')
% save(outVarFullName ,'outVar')



%     % Single out the M5.1 event from SCSN with a lot of records
%     is_scsn  = (strcmp(trList.dataSetName,'scsn')|strcmp(trList.dataSetName,'scsnPx'));
%     has_m5p1 = (trList.eq.m>5.09 & trList.eq.m<5.11);
%     idx      = find(is_scsn & has_m5p1);
%     tmpList  = trList.selectSubList(idx);
%     
%     idx_woM5p1 = setdiff(1:numel(trList.eq.m),idx);
%     testList   = trList.selectSubList(idx_woM5p1); 
%     trList     = testList;
%     
%     % Plot wforms
%     ntrMax = 10;    % From each mBin randomly select <ntrMax> traces, 
%     nruns  = 5;     % repeat <nruns> times
%     [hf2]  = plot_wforms_for_each_range(trList,mRanges,rRange,tRange,ntrMax,nruns,fig);
%     
%     % Plot data
%     trList   = testList;
%     idxPlot  = find(trList.eq.m>min(mRanges(:)) & trList.dist.hyp<=max(rRange(:)));
%     plotList = trList.selectSubList(idxPlot);
%     summary_trList(plotList,[],1,1);
% 
%     % Plot peak amplitude movie
%     plot_peakAmps_movie(trList,mRanges,rRange,tRange,band)
%     plot_peakAmps_movie(zList,mRanges,rRange,tRange,band)
% 
%     % Compute and plot tsep matrix
%     fig.print = false;
%     plevel  = 0.05;
%     tRange  = [0, 20];
%     rRange  = [0 25];
%     mRanges = [6.5 7.7; 6 6.5; 5.5 6; 5 5.5; 4.5 5; 4 4.5; 3.5 4; 3 3.5; 2.5 3];
%     [hf2]   = get_tsep_matrix(nsList,'Pd',mRanges,rRange,plevel,fig);
%     
%     
% 
%     %% Plot single ground motion curve
%     i = 0;
%     i = i+1;
%     singleIdx  = nchi;
%     idx_chichi = find(zList.eq.m>7 & zList.dist.hyp<20 & zList.snr>=snrMin);
%     chiList    = zList.selectSubList(idx_chichi);
%     [hf]       = plot_single_PGx_curve(chiList,singleIdx,[-5 30]);
%     
%     % South Napa records
%     [val,sortIdx] = sort(zList.dist.hyp);
%     zList.sortList(sortIdx)
%     [hf] = plot_single_PGx_curve(zList,1,[-5 30]);
%     
%     idx_m7 = find(zList.eq.m>7 & zList.dist.hyp<50);
%     m7List = zList.selectSubList(idx_m7);
%     % Good examples in m7List: 76, 80, 81 (r=50), 82(r=35)
%     [hf]   = plot_single_PGx_curve(m7List,76,[-5 50]);
%     
%     nchi = numel(idx_chichi);
%     for ichi = 1:nchi
%     
%         [S,meta] = read_any_trace(chiList.fullName{ichi},chiList,1);
%         ppxIdx   = chiList.ppxIdx(ichi);
%         h=plot_wform(S.vel,meta.t,meta.t(ppxIdx),[],1);
%         chiList.snr(ichi)
%         pause
%     end
% 
%     
%     
%     %% Plot tauC and tauPmax 
%     fprintf(1,'\nFor now, use i33-outfiles for plotting tauPmax and tauC.\n Recompute paramInf and re-estimate them in i34!\n\n')
%     load '/scratch/memeier/fbout/i33/causal_0p5/SSR_marginalPost_noWt/zList_paramErrors.mat'
%     
%     % Full List
%     tmpList = out.zList;
%     tmpEqs  = compile_eqs_structure(tmpList);
%     [hfig]  = plot_tauPC(tmpList,tmpEqs);
% 
%     %mRanges = [6.5 7.7; 6 6.5; 5.5 6; 5 5.5; 4.5 5; 4 4.5; 3 4];
%     mRanges = [6.5 7.7; 6 6.5; 5.5 6; 5 5.5; 4.5 5; 4 4.5];
%     rRange  = [0 25];
%     idx_tmp = find( tmpList.eq.m>=min(mRanges(:)) & tmpList.eq.m<max(mRanges(:)) & tmpList.dist.hyp>=rRange(1) & tmpList.dist.hyp<rRange(2) );
%     nsList  = tmpList.selectSubList(idx_tmp);
%     nsEqs   = compile_eqs_structure(nsList);
%     
%     [hfig]  = plot_tauPC(nsList,nsEqs);
%     %trList     = hList;                 % hList; zList;
%     %trList     = zList;
%     %trList     = kList;
%     %trList     = kikList;
%     
%     %     [hf]   = cf_peak_amps(trList,'Pa',[tmin tmax],pctl,yscale,o_printFig);
%     %     for iband = 1:nbands
%     %         [hf] = cf_peak_amps(trList,iband,[tmin tmax],pctl,yscale,o_printFig);
%     %     end
%     
%     % Add tdpa curves from individual events to plot
%         % Peak amps
%         
% %    idx_chichi=find(zList.eq.m>7.6 & zList.eq.m<7.7 & zList.dist.hyp<25 & zList.snr>=snrMin);
%     %idx_chichi = find(zList.eq.m>7 & zList.dist.hyp<25 & zList.snr>=snrMin);
%     %chiList    = hList.selectSubList(idx_chichi);
%     idx_chichi = find(hList.eq.m>7 & hList.dist.hyp<25 & hList.snr>=snrMin);
%     chiList    = hList.selectSubList(idx_chichi);
%     amps       = chiList.fb.amax;
%     pd         = chiList.pd;
%     
%     % Fill in zeros until all time series are as long as longest one
%     nmax = max(cellfun(@(x) numel(x), pd));
%     ntr  = numel(pd);
%     Pd   = zeros(ntr,nmax);
%     
%     for itr = 1:ntr
%         na        = length(pd{itr});
%         Pd(itr,:) = [pd{itr}, repmat(pd{itr}(end),1,nmax-na)];
%     end
%     
%     xlms=get(gca,'xlim');
%     t = linspace(xlms(1),xlms(2),nmax);
%     
%     nlines = 4 + size(mRanges,1)-1;
%     figure(hf); subplot(nlines,6,[1:4,7:10,13:16,19:22]); hold on;
%     plot(t,log10(Pd'),'b','lineWidth',2)
%     
% 
%     [signal{ir},noise{ir},idx{ir}] = ...
%         get_time_dependent_peak_amps(trList,mRanges(ir,:),rRange,band,pctl,fig.yscale);
% 
%     
%     %% Do m5+ have lower absolute noise amps than m6+? YES! goddamit! why?
%     [nAcc5,nVel5,nDsp5] = get_absolute_noise_amps(trList,[5 6]);
%     [nAcc6,nVel6,nDsp6] = get_absolute_noise_amps(trList,[6 7]);
%     
%     noisePct5 = prctile(nAcc5,[5 50 95]);
%     noisePct6 = prctile(nAcc6,[5 50 95]);
%     
%     figure(192); clf; hold on; grid on;
%     set(gca,'yscale','log','fontSize',ftSize)
%     p1 = plot(noisePct5','color','r','lineWidth',2);
%     p2 = plot(noisePct6','color','k','lineWidth',2);
%     
%     l1 = legend([p1(1);p2(1)],'m5+ (5^{th}, 50^{th} & 95^{th} percentiles)','m6+ (5^{th}, 50^{th} & 95^{th} percentiles)');
%     set(l1,'fontSize',ftSize)
%     xTix = fliplr(get(gca,'xtick'))*snpLength;
%     set(gca,'xtickLabel',xTix,'fontSize',ftSize)
%     xlabel('Time before p-pick [sec]','fontSize',ftSize)
%     ylabel('Absolute GM levels [SI-units]','fontSize',ftSize)
%     title(sprintf('Absolute noise levels for all traces with SNR>%i',snrMin),'fontSize',ftSize)
%     
%     if o.printNoise
%         figPath     = sprintf('~/programs/seismo/fig/i%i/snr/new/',iN);
%         figFullName = sprintf('%snoiseAmps_snrMin%i',figPath,snrMin);
%         set(gcf,'PaperPositionMode','auto')
%         print('-depsc2',[figFullName,'.eps'])
%         print('-dpng',[figFullName,'.png'])
%     end


% if o.correctPickDelay
% 
% 
%     %clf; plot(tmp(:,1),log10(tmp(:,2)),'xk')
%     %set(gca,'xlim',[-1 15],'ylim',[-7 -2])
%     
%     fprintf(1,'Correcting for pick delay ... ')
%     idx    = find(zList.dist.hyp<25);
% 	nsList = zList.selectSubList(idx);
%     %nsList = zList;
%     nz     = numel(nsList.eq.m);
%     
%     tmp    = cell2mat(nsList.var2);
%     deltas = tmp(:,1);
%     nAmps  = tmp(:,2);
%     
%     lgc = (deltas~=99999);
%     idx = find(lgc);
%     nz  = sum(lgc);
%     
%     nsList = nsList.selectSubList(idx);
%     deltas = deltas(lgc);
%     nAmps  = nAmps(lgc);
%     pd_bak = nsList.pd;
%     
%     
%     for iz = 1:nz
%         
%         ds            = deltas(iz);
%         pd            = nsList.pd{iz};
%         pd_shifted    = [zeros(1,ds), pd];
%         nsList.pd{iz} = pd_shifted;
%         
%         %         if ds==99999
%         %             fprinf(1,'recompute delay with shorter picker windows. do it. now\n')
%         %             fprinf(1,'set ds for manual picks to zero?\n')
%         %             pause
%         %             
%         %             
%         %         elseif ds~=0
%         %             pd            = nsList.pd{iz};
%         %             pd_shifted    = [zeros(1,ds), pd];
%         %             nsList.pd{iz} = pd_shifted;
%         %         end
%     end
%     fprintf(1,'done.\n')
% end
% 
% 
% % Recompute pick shift
% o.neverTrue = false;
% if o.neverTrue
%     n = 2.5;
%     k = 1;
%     
%     for iz = 1:nz
%         % Read waveform
%         % Process trace
%         % Compute shift
%         % Reprocess waveform with new pick
%         % Compute pd with new pick
%     end
% end



% tPrimeSignal = .08;       % Measure signal amp at tPrimeSignal
% tPrimeNoise  = 1;       % Measure noise amp at tPrimeNoise
% nsnoise = cellfun(@(x) numel(x), zList.accNoise);
% sr      = zList.sRate;
% m       = zList.eq.m;
% iPrimeNoise                      = nsnoise - tPrimeNoise*sr;    % Noise sample <tPrimeNoise> sec before pick
% iPrimeNoise(iPrimeNoise<1)       = 50;                          % If that is negative, set it to 50
% iPrimeNoise(iPrimeNoise>nsnoise) = 1;                           % If 50 is more than there are noise samples, set it to 1
% iPrimeSignal = ceil(tPrimeSignal*sr);
% nz    = numel(zList.eq.m);
% naAmp = zeros(nz,1);
% saAmp = zeros(nz,1);
% for i = 1:nz
% naAmp(i) = zList.accNoise{i}(iPrimeNoise(i));
% saAmp(i) = zList.pa{i}(iPrimeSignal(i));
% end
% clf; hold on; 
% sc = scatter(naAmp,saAmp,12,m,'o','filled');
% set(gca,'xscale','log','yscale','log')
% grid on
% colorbar
% set(gca,'xlim',[1e-6 1],'ylim',[1e-7 1],'fontSize',ftSize)
% xlabel(sprintf('Noise amps at %5.2fsec before pick [m/s/s]',tPrimeNoise),'fontSize',ftSize)
% ylabel(sprintf('Signal amps at %5.2fsec after pick [m/s/s]',tPrimeSignal),'fontSize',ftSize)
% 
% set(gcf,'PaperPositionMode','auto')
% print('-depsc','../../fig/i35/var/noiseAmps_vs_signalAmps_4.eps')
