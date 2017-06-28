function [trList] = recompute_Pd_with_PDC(trList,params,opt)

global iN fOrder snpLength

hundreds  = linspace(1e2,1e5,1e3);
thousands = linspace(1e3,1e6,1e3);

fprintf(1,'8ung: some hardcoded parameters. check.\n')

expVect = params.expVect;
k       = params.k;
nrnd    = params.nrnd;
fMode   = opt.fMode;

ntmp    = numel(trList.m);
nexp    = numel(expVect);

trList.var3 = cell        (ntmp,nexp);      % Use for PD matrix
trList.var4 = single(zeros(ntmp,nexp));     % Use for DS matrix


% PICKER (see SBPx.m for meaning of parameters)
px           = load_picker_settings('production_short');
fLow_prefilt = px.Param.fLow_prefilt;
fLow_px      = px.Param.fLow_px;
fUp_px       = px.Param.fUp_px;
ntap         = px.Param.ntap;

tnoise   = 4;              % Measure peak amps starting <tnoise> sec before p-pick
tmax     = 20;             % Stop measuring peak amps after <tmax> sec
tmax_bak = tmax;
nsnpmin  = 500;

% Compose file names for saving modified trList and additional
% information
outDirName     = sprintf('/scratch/memeier/fbout/i%i/zLists/',iN);
exponentString = strrep(sprintf('%3.1f_',expVect),'.','p');
outName        = sprintf('ntr%i_n%s',ntmp,exponentString(1:end-1));
outDirFullName = sprintf('%s%s/',outDirName,outName);
outFullName    = sprintf('%s%s',outDirFullName,'trList.mat');

if ~exist(outDirFullName,'dir'); fprintf(1,sprintf('Created directory %s\n',outDirFullName))
    unix(['mkdir -p ',outDirFullName]);
else                             fprintf(1,'\nWARNING: out-directory for this modified Pd-list already exists.\n')
    fprintf(1,'         If you proceed, all existing file will be overwritten.\n')
    fprintf(1,'         Think twice.\n'); pause
end


% Loop over all traces
tic
for itr=1:ntmp
    
    %if (ismember(itr,thousands) |itr==1); fprintf(1,sprintf('\n%i/%i\n',itr,ntmp)); end
    if (ismember(itr,hundreds) |itr==1); fprintf(1,sprintf('\n%i/%i\n',itr,ntmp)); end
    
    % Read waveform and meta info
    [S,meta]  = read_any_trace(trList.fullName{itr},trList,opt);
    m         = trList.m(itr);
    t         = meta.t;
    sr        = meta.sr;
    ppxIdx    = trList.ppxIdx(itr);          % Read original, uncorrected pick
    instrCode = trList.instrCode(itr);
    isSM      = (strcmp(instrCode,'L') ||strcmp(instrCode,'N') ||strcmp(instrCode,'G'));
    
    % Prepare vertical trace for picking
    smpx    = S.raw - mean(S.raw(1:2*ntap));                              % Remove mean
    stappx  = taper(smpx,sr,ntap);                                        % Taper
    if ~isSM; [sh_px]   = bworth(stappx,sr,fLow_prefilt ,'high',2,'causal');          % Pre-filter (high-pass)
              [apx,~,~] = bb2accVelDsp(sh_px,sr,ppxIdx,'allWform',fLow_prefilt,2,fMode);
              abpx      = bworth(apx   ,sr,[fLow_px fUp_px],'band',2,'causal');
    else      abpx      = bworth(stappx,sr,[fLow_px fUp_px],'band',2,'causal');
    end
    
    % Compute median pick delay
    metapx.t      = t;
    metapx.ppxIdx = ppxIdx;
    metapx.sr     = sr;
    
    for iexp = 1:nexp
        
        % Exponent
        n = expVect(iexp);
        if n==0; ds=0;
        else     [deltasVect ,maxNoise ] = get_median_pick_delay(abpx,metapx,n,k,nrnd,px,opt);
            ds                           = round(median(deltasVect (deltasVect ~=99999)));
        end
        
        % If pick delay is reasonable reprocess waveform with new/corrected pick
        trList.var4(itr,iexp) = ds;
        if ~isnan(ds);
            
            % Corrected pick
            newppxIdx = ppxIdx-ds;
            
            % Remove early mean, taper and pre-filter waveform (2-pole high pass) ...
            sm   = S.raw - mean(S.raw(1:newppxIdx));
            stap = taper(sm,sr,ntap);
            sb   = bworth(stap,sr,[fLow_prefilt,30],'band',2,fMode);
            if isSM; [acc2,vel2,dsp2] = sm2accVelDsp(sb,sr,newppxIdx,'afterPx',fLow_prefilt,2,fMode);
            else     [acc2,vel2,dsp2] = bb2accVelDsp(sb,sr,newppxIdx,'afterPx',fLow_prefilt,2,fMode);
            end
            
            % Recompute pd-curve
            if m>7; tmax = 60; end
            [pd,pnd,pdIdx,~,~]    = measure_tdpa(dsp2,newppxIdx,sr,tnoise,tmax,snpLength,[]);
            %[pv,pnv,pvIdx,~,~]    = measure_tdpa(vel2,newppxIdx,sr,tnoise,tmax,snpLength,[]);
            %[pa,pna,paIdx,~,~]    = measure_tdpa(acc2,newppxIdx,sr,tnoise,tmax,snpLength,[]);
            pd                    = trim_amax(pd,nsnpmin);
            trList.var3{itr,iexp} = single(pd);
            if m>7; tmax = tmax_bak; end    % Reset to original value
            %tppx = t(newppxIdx);plot_wform_and_tdpa(dsp2,t,tppx,[],[tppx-1.1*tnoise 1.1*tppx+tmax],pd,snpLength,1115)
        end
    end
    
    % Save temporary trList
    if ismember(itr,thousands) fprintf(1,sprintf('\n%i/%i ... saving temporary outfile ...',itr,ntmp))
        outTmpFullName = strrep(outFullName,'.mat',sprintf('_%i.mat',itr));
        save(outTmpFullName,'trList'); fprintf(1,'done. '); toc
    end
end


if opt.saveOut
    
    % Save final version of trList
    outListFullName = sprintf('%strList.mat',outDirFullName);
    save(outListFullName,'trList')
    %save('/scratch/memeier/fbout/i36/zLists/ntr914_n0p0_1p0_2p0_3p0_4p0/nssubList_914.mat','trList')
    
    % Save complementary information
    outVar.expVect          = expVect;
    outVar.dataSetBaseNames = dataSetBaseNames;
    outVar.useBlackList     = o.useBlackList;
    outVar.usePxList        = o.usePxList;
    outVar.configFileName   = configFileName;
    outVar.saveTimeStamp    = clock;
    outVar.snpLength        = snpLength;
    outVar.nrnd             = nrnd;
    outVarFullName = sprintf('%s%s',outDirFullName,'var.mat');
    save(outVarFullName ,'outVar')
end
toc