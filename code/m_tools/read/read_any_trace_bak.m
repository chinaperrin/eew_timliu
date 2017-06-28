function [S,meta] = read_any_trace(traceFullName,trList,o)

% Reads gain-corrected Cali and un-corrected Japan waveforms, applies
% baseline correction, tapering, high-pass pre-filtering and - if the
% waveform is an accelergram - integrates to velocity. Returns acceleration
% velocity and displacement. Assumes that there is a p-pick saved in trList.

% For Japan:    gain correction is applied
% For NGA:      dito
% For SoCal:    expects input trace to be already gain corrected


% OPTIONS
% if o.process = 1           mean is removed, tapered and highpass
%                            filtered, integrated and differentiated to
%                            acc, vel & dsp
% if o.process = 0           raw waveform is returned
%
% if o.intMode = 'allWform'  trace is integrated from first sample on
% if o.intMod  = 'afterPx'   integration starts at pick-index
%
%
% ISSUE: not sure if I should highpass integrated displacement traces. From
%        the examples I have seen it sometimes helps, sometimes makes it
%        worse ...
%
%
% WAVEFORM VARIABLE NAMES
% First letter:         Nature of raw traces 
%                       a       acceleration (SM record)
%                       v       velocity (BB record)
%
% Following letters:    Applied processing steps
%                       h,l,b   applied filters (high-, low-, bandpass)
%                       i       integration
%                       d       differentiation
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .


addpath(genpath('../../../matlab/sac/')) 
addpath(genpath('../m_filter/'))

if (iscell(traceFullName)); error('Give traceFullName as a string, not as a cell'); end


% Download waveform from bigstar if necessary
if ( (~exist(traceFullName,'file')) && (o.scp_wforms) )
    scp_wform(traceFullName)
end
    

% Filter settings for high-pass pre-filtering
fOrder       = 4;
fMode        = 'causal';
fLow_prefilt = 0.075;
ntap         = 100;
%fprintf(1,sprintf('Using %s filter of order %i... appropriate?\n', fMode, fOrder))

trIdx     = find(strcmp(traceFullName,trList.fullName));    % Find trace in trList
ptIdx     = regexp(traceFullName,'\.');                     % Find out what sort of waveform it is
extension = traceFullName(ptIdx(end):end);

jp_extensions = {'.EW','.EW1','.EW2','.NS','.NS1','.NS2','.UD','.UD1','.UD2'};
 

% JAPAN - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ismember(extension,jp_extensions)
    
    % Read waveform
    [sraw,meta] = read_ascii_wform_jp(traceFullName,17,1);
    ppxIdx      = trList.ppxIdx(trIdx);
    sr          = meta.sr;

    % Save meta info
    meta.ppxIdx  = ppxIdx;
    meta.t       = 0:1/sr:(numel(sraw)-1)/sr;
    meta. origin = 'jp';
    
    % Apply gain correction
    slashIdx     = find(meta.scaleFactor=='/');
    scaleFactor1 = str2double(meta.scaleFactor(slashIdx+1:end));
    bracketIdx   = find(meta.scaleFactor=='(');
    scaleFactor2 = str2double(meta.scaleFactor(1:bracketIdx-1));
    S.raw        = sraw*scaleFactor2/scaleFactor1*0.01;                    % From [counts] to [gal] to [ms^-2]
    
    if o.process
        
        % If pick is either empty or = 0
        if isempty(ppxIdx); ppxIdx = 0; end
        if ppxIdx==0
            fprintf(1,'No p-pick available. Using mean of t = [1:2*ntap] for mean removal\n')
            ppxIdx = 2*ntap;
        end
        
        sm     = S.raw - mean(S.raw(1:ppxIdx));                           % Remove early mean from gain corrected wform...
        stap   = taper(sm,sr,ntap);                                       % Taper waveform ...
        ah     = bworth(stap,sr,fLow_prefilt,'high',fOrder,fMode);        % Prefilter ...
        ahi    = integrate_wform(ah,sr,ppxIdx,o.intMode);                 % Integrate to vel (from [m/s^2] to [m/s])
        ahih   = bworth(ahi,sr,fLow_prefilt,'high',fOrder,fMode);         % Highpass
        ahihi  = integrate_wform(ahih,sr,ppxIdx,o.intMode);               % Integrate to dsp (from [m/s] to [m])
        ahihih = bworth(ahihi,sr,fLow_prefilt,'high',fOrder,fMode);       % Highpass

        S.acc  = ah;
        S.vel  = ahih;
        S.dsp  = ahihih;
        S.dsp2 = ahihi;
    end
    
    
    
% CALIFORNIA  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif strcmp(extension,'.sac')
    
    [out] = read_sac_trace3(traceFullName);
    sr    = out.sr;
    t     = out.t;
    S.raw = out.sraw/100;                                                  % From [cm/s] to [m/s] or
    
    if max(S.raw)>10000; 
        fprintf(1,'Maximum amplitude is >1e4. Are you sure it is a gain corrected trace?\n')
        pause
    end
    
    meta.t      = t;
    meta.sr     = sr;
    meta.origin = 'cali';
    
    % If there is no issue with the trace ...
    if (sr~=0)
        
        tppx           = trList.tppx(trIdx);
        [~,ppxIdx]     = min(abs(tppx - t));                               % The index of the pick
        meta.ppxIdx    = ppxIdx;
        meta.skipTrace = false;
    
        if o.process
            if tppx==0
                fprintf(1,'No p-pick available. Randomly assigning ppxIdx = 2*ntap...\n Integration only after ppxIdx ...\n')
                ppxIdx = 2*ntap;
            end
            
            sm   = S.raw - mean(S.raw(1:ppxIdx));                          % Remove early mean ...
            stap = taper(sm,sr,ntap);                                      % Taper waveform ...
            [sh] = bworth(stap,sr,fLow_prefilt,'high',fOrder,fMode);       % Prefilter
            
            % Seismo- or accelerogram
            traceSpecs = textscan(traceFullName,'%s','delimiter','/');
            traceName  = traceSpecs{1}{end};
            traceInfo  = textscan(traceName,'%s','delimiter','.');
            chan       = traceInfo{1}{4};
            
            % Strong motion records
            if ( strcmp(chan(2),'L') | strcmp(chan(2),'N') | strcmp(chan(2),'G') )

                % Could also use function sm2accVelDsp.m
                ah      = sh;                                                % Highpassed raw signal = acceleration
                ahi     = integrate_wform(ah,sr,ppxIdx,o.intMode);           % Integrate to vel (from [m/s^2] to [m/s])
                ahih    = bworth(ahi,sr,fLow_prefilt,'high',fOrder,fMode);   % Highpass
                ahihi   = integrate_wform(ahih,sr,ppxIdx,o.intMode);         % Integrate to dsp (from [m/s] to [m])
                ahihih  = bworth(ahihi,sr,fLow_prefilt,'high',fOrder,fMode); % Highpass
                
                S.acc  = ah;                                                  % Assign to S
                S.vel  = ahih;      
                S.dsp  = ahihih;  
                S.dsp2 = ahihi;
                
                
            % Broad band records
            elseif (strcmp(chan(2),'H'))
                
                % Could also use function bb2accVelDsp.m
                vh    = sh;                                                 % Highpassed raw signal = velocity
                vhd   = diff(vh)*sr;                                        % Differentiate to acc (from [m/s] to [m/s^2])
                vhd   = [vhd; vhd(end)];                                    % Duplicate last sample to have equal length
                vhi   = integrate_wform(vh,sr,ppxIdx,o.intMode);            % Integrate to dsp (from [m/s] to [m])
                vhih  = bworth(vhi,sr,fLow_prefilt,'high',fOrder,fMode);    % Highpass
                
                S.acc  = vhd;                                                % Assign to S
                S.vel  = vh;
                S.dsp  = vhih;
                S.dsp2 = vhi;    
            end
        end
    end
    
    
% NGA West1   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif strcmp(extension,'.AT2')
    
        [sraw,meta] = read_ascii_wform_nga(traceFullName,1);
            
        % From [g] to [ms^-2] (BIPM: 1g = 9.80665 ms^-2)
        S.raw = sraw*9.80665;
    
        % Meta info
        sr           = meta.sr;
        ppxIdx       = trList.ppxIdx(trIdx);                               % The index of the pick
        meta.ppxIdx  = ppxIdx;
        meta.t       = 0:1/sr:(numel(sraw)-1)/sr;
        meta. origin = 'ngawest';

        if (ppxIdx~=0); sm = S.raw-mean(S.raw(1:ppxIdx));                  % Remove early mean ...
        else            sm = S.raw-mean(S.raw(1:2*ntap));
        end
 
        if o.process 

            stap   = taper(sm,sr,ntap);                                   % Taper waveform ...
            ah     = bworth(stap,sr,fLow_prefilt,'high',fOrder,fMode);     % Prefilter ...
            ahi    = integrate_wform(ah,sr,ppxIdx,o.intMode);              % Integrate to vel (from [m/s^2] to [m/s])
            ahih   = bworth(ahi,sr,fLow_prefilt,'high',fOrder,fMode);      % Highpass
            ahihi  = integrate_wform(ahih,sr,ppxIdx,o.intMode);            % Integrate to dsp (from [m/s] to [m])
            ahihih = bworth(ahihi,sr,fLow_prefilt,'high',fOrder,fMode);    % Highpass
            
            S.acc  = ah;
            S.vel  = ahih;
            S.dsp  = ahihih;
            S.dsp2 = ahihi;
        end

        
    
% TURKEY      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
elseif strcmp(extension,'.txt')

%     [sraw,meta] = read_ascii_wform_tk(traceFullName,1,1);
% 
%     fprintf(1,'8ung: current way of converting units is approximate. check.\n')
%     pause
% 
%     if (strcmp(ornt,'E'))           % Choose trace and convert from [gal] to [m/s^2]
%         sraw = sraw.ew/100;    
%     elseif (strcmp(ornt,'N'))
%         sraw = sraw.ns/100;
%     elseif (strcmp(ornt,'Z'))
%         sraw = sraw.ud/100;
%     end    
%     
%     ppxIdx   = trList.ppxIdx(trIdx(1));
%     sr      = meta.sr;    
%     
%     if (ppxIdx>0)
%         sm  = sraw - mean(sraw(1:ppxIdx));                                  % Remove early mean ...
%     else
%         fprintf(1,'No pick found\n')
%         sm  = sraw - mean(sraw(1:3*ntap)); 
%     end
%     stap    = taper(sm,sr,ntap);                                           % Taper waveform ...
%     [s,~,~] = butter_pass_tdomain_f(stap,fLow_prefilt,999,sr,7,npad,0);    % Prefilter ...                                       % Integrate to velocity  ...
%     
%     S.acc   = s;
%     [S.vel] = acc2vel(s,sr,2);      % Integrate to velocity     (from [m/s^2] to [m/s])
%     [S.dsp] = acc2vel(S.vel,sr,2);  % Integrate to displacement (from [m/s]   to [m])
%     
%     meta.t      = 0:1/sr:(numel(s)-1)/sr;
%     meta.sr     = sr;
%     meta.ppxIdx  = ppxIdx;
%     meta.origin = 'turkey';

else
    fprintf(1,'EXTENSION NOT FOUND, WHAT NOW?\n')
    pause
end