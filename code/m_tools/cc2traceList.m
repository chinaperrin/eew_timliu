clear all
% Reads in Carlo Cauzzi's strong motion data set and turns it into a
% traceList compatible with filterBank_xxx.m
% menandrin@gmail.com, 130424



% ISSUES
% - do horizontal records always have same number of header lines?
% - same nr of columns?
% - do the multi-column ascii waveforms sometimes contain empty entries in
%   last columns?
% - are they all in cm/s^2? prob. m/s^2
% - normal to get low-freq wiggles upon integration
% - already some processing done? mean removed?
% - is real time integration feasible?
% - how does cc handle wrong/inaccurate picks?
% - can cc explain the low freq transients?

o_saveTraceList   = 0;

wformDir          = '/scratch/memeier/VS/data/wform/cc/wforms/raw/';
TraceListFileName = strcat([wformDir,'trList.mat']);
matFileName       = '/scratch/memeier/VS/data/wform/cc/meta/CF08_databank_V_added.mat';




fprintf('\n\n\n\n\n\n###########################################################\n')
fprintf('                    cc2traceList\n')
fprintf('###########################################################\n\n')

% Load meta-data file (complied on Yavor's computeron xXXX)
load(matFileName); clear num txt
header = raw(1,:);
raw    = raw(2:end,:);

idx    = [(888:961)';(988:1005)';(1040:1080)']
data   = raw(idx,:);

noNewEq     = cell2mat(cellfun(@(x) isnan(x(:,1)), raw(:,1),'uniformOutput',false));
lastEqIdx   = find(noNewEq==0,1,'last'); % --> entries which are traces, not new earthquakes
hasVertical = cell2mat(raw(:,24));
%noProc     = strcmp(correctionNote,'pre-event BC');

nlines    = size(raw,1);
traceList = cell(nlines,26);
traceIdx  = 0;
ieq       = 0;


for iline = 1:nlines
    
    if (~noNewEq(iline))
        
        fprintf(1,['\nEarthquake: ',raw{iline,1}])
        ieq = ieq + 1;
        
        eqNameString = raw{iline,1};
        ulineIdx     = regexp(eqNameString,'_');
        dateStr      = eqNameString(1:ulineIdx(3)-1);
        eqDate       = strrep(dateStr, '_', '/');
        t0           = eqNameString(ulineIdx(3)+1:end);
       
        Mw     = raw{iline,6};
        eqMech = raw{iline,3};
        eqLat  = raw{iline,5};
        eqLon  = raw{iline,4};
        eqZ    = raw{iline,7};
        
        
        % Count stations which have recorded this event
        if (iline < lastEqIdx)
            tmp = noNewEq(iline+1:end);
            nst = find(tmp==0,1,'first');
        else
            nst = nlines-iline+1;
        end
        fprintf(1,[' has recordings from ',num2str(nst),' stations\n'])
        

        % Extract info from each station on which there are records of this earthquake
        for ist = iline:iline+nst-1
            fprintf(1,['  station ',num2str(raw{ist,8}),'\n'])
            
            if (hasVertical(ist)==1)
                ntr = 3;
            else
                ntr = 2;
            end
            
            % Go through all traces on this station
            for itr = 1:ntr
                
                % ... add trace to list
                traceIdx = traceIdx + 1;
                
                if (itr == 1)
                    
                    % First horizontal component
                    wformName   = strcat(wformDir,raw{ist,16});
                    nlHdr       = raw{ist,11};
                    nCol        = raw{ist,12};
                    corrComment = raw{ist,21};
                    
                elseif (itr == 2)
                    
                    % Second horizontal component
                    wformName   = strcat(wformDir,raw{ist,17});
                	nlHdr       = raw{ist,11};
                    nCol        = raw{ist,12};
                    corrComment = raw{ist,21};
                    
                elseif (itr == 3)
                    
                    % Vertical component
                    wformName   = strcat(wformDir,raw{ist,25});
                    nlHdr       = raw{ist,27};
                    nCol        = raw{ist,28};
                    corrComment = raw{ist,26};
                    
                end
                
                fprintf(1,['        ',wformName,'\n'])
                
                % Write to traceList
                traceList{traceIdx,1}  = wformName;
                traceList{traceIdx,2}  = t0;
                traceList{traceIdx,4}  = raw{ist,18};     % epicentral distance
                traceList{traceIdx,7}  = eqLat;
                traceList{traceIdx,8}  = eqLon;
                traceList{traceIdx,9}  = eqZ;
                traceList{traceIdx,10} = eqDate;
                
                traceList{traceIdx,11} = 'cc';              % data set name
                traceList{traceIdx,12} = raw{ist,8};        % staCode
                traceList{traceIdx,13} = raw{ist,9};        % staName
                traceList{traceIdx,14} = raw{ist,22};       % stLat
                traceList{traceIdx,15} = raw{ist,23};       % stLon
                traceList{traceIdx,16} = raw{ist,10};       % sr
                traceList{traceIdx,17} = Mw;
                traceList{traceIdx,18} = eqMech;
                traceList{traceIdx,19} = nlHdr;             % Number of header lines
                traceList{traceIdx,20} = nCol;              % Number of columns in wform file
                traceList{traceIdx,21} = raw{ist,13};       % attributed ground category
                traceList{traceIdx,22} = raw{ist,14};       % Instr
                traceList{traceIdx,23} = raw{ist,15};       % rawDataSrc
                traceList{traceIdx,24} = corrComment;       % comment on wform correction
                traceList{traceIdx,25} = raw{ist,19};       % ipoDist   
                traceList{traceIdx,26} = raw{ist,20};       % faultDist 
       
            end
            

        end
        1+1;
    end
end
fprintf(1,[' \n\n ',num2str(traceIdx), ' traces from ',num2str(ieq),' earthquakes have been extracted. Thanks Carlo, great job.\n'])


% Save traceList
if (o_saveTraceList == 1)
    fprintf(1,'Saving TraceList ...')
    save(TraceListFileName,'traceList');
    fprintf(1,' DONE\n')
end