function [stacks,counts] = stack_wforms_around_px(trList,ns_prePx,ns_postPx,opts)

% Stacks absolute values of waveform traces around their pick. Downsamples
% traces with sRate=200 to sRate=100.
%
% ns_prePx  = no. of samples before pick
% ns_postPx = no. of samples after  pick
%
% menandrin@gmail.com, 150816


% Current limitations
% - uses acceleration traces (hard-coded)
% - skips traces if they are shorter than requested interval (--> argin)
% - skips all traces that have sRate ~=100 or 200

fprintf(1,'This is an old function. There should be a better one by now that
uses interpolation to get a common time sampling for all traces. Consider using
that one instead\n')
pause

ntr = numel(trList.m);

ns_stack  = ns_prePx + 1 + ns_postPx;
accStack  = zeros(ntr,ns_stack);
velStack  = zeros(ntr,ns_stack);
dspStack  = zeros(ntr,ns_stack);
idx_stack = -ns_prePx:1:ns_postPx;

flg_wrongSR       = false(ntr,1);
flg_shortSignal   = false(ntr,1);
flg_wformNotFound = false(ntr,1);

tens     = linspace(1e1,1e4,1e3);
hundreds = linspace(1e2,1e5,1e3);

if opts.verbose; figure(1); clf; hold on; grid on; box on; end



%for itr = 1:ntr
parfor itr = 1:ntr
    
    fullName      = trList.fullName{itr};
    flg_skipTrace = false;
    absAcc        = 0;
    absVel        = 0;
    absDsp        = 0;
    
    if opts.verbose
        [~,recordName,~] = splitFullName(fullName);
        fprintf(1,sprintf('\n%i/%i: \trecord %s\n',itr,ntr,recordName))
        fprintf(1,['\t\t',num2str(trList.m(itr),'%3.1f'),'M @ ',num2str(trList.hypDist(itr),'%3.0f'),'km\n'])
        fprintf(1,sprintf('\t\tSNR: %12.0f\n ',trList.snr(itr)))
    else
        if ismember(itr,tens); fprintf(1,sprintf('%i/%i\n',itr,ntr)); end
    end
    
    % Copy wform from bigstar if necessary
    if ( (~exist(trList.fullName{itr},'file')) && (opts.scp_wforms) );
        fprintf(1,'scp-ing wform\n')
        scp_wform(fullName);
    elseif ( (~exist(trList.fullName{itr},'file')) && ~(opts.scp_wforms) );
        fprintf(1,'Wform not available and opts.scp_wforms set to zero; skipping this trace\n')
        flg_skipTrace          = true;
        flg_wformNotFound(itr) = true;
    end
    
    % Read trace
    [s,meta] = read_any_trace(trList.fullName{itr},trList,opts);
    t        = meta.t;
    ppxIdx   = meta.ppxIdx;
    tppx     = t(ppxIdx);
    sr       = meta.sr;
    ns       = numel(s.raw);
    
    if sr==100        
        % Use signal as is
        absAcc = abs(s.acc);
        absVel = abs(s.vel);
        absDsp = abs(s.dsp);
        %time   = t;
        
    elseif sr==200
        % Downsample
        if (mod(ppxIdx,2)==0);                          % Even ppxIdx
            dsIdx     = 2:2:ns;
            ppxIdx = ppxIdx/2;
            %fprintf(1,' Even, ')
        else                                            % Odd  ppxIdx
            dsIdx     = 1:2:ns;
            ppxIdx = ppxIdx/2+.5;
            %fprintf(1,' Odd,  ')
        end
        absAcc = abs(s.acc(dsIdx));
        absVel = abs(s.vel(dsIdx));
        absDsp = abs(s.dsp(dsIdx));
        %time   = t        (dsIdx);
        
        if opts.verbose; dt = time(ppxIdx)-tppx;
                         fprintf(1,sprintf('\tdt = %6.3f\n',dt))
        end
        
    else
        fprintf(1,'SR not equal 100 or 200; skipping this trace\n')
        flg_skipTrace    = true;
        flg_wrongSR(itr) = true;
    end
    
   
    
    % SAVE VALUES AROUND PPX IN STACK
    if ~flg_skipTrace
        
        ns   = numel(absAcc);
        sIdx = ppxIdx - ns_prePx;
        eIdx = ppxIdx + ns_postPx;
        
        if ( sIdx>0 & eIdx<=ns )
            
            % Add up signals
            accStack(itr,:) = absAcc(sIdx:eIdx);
            velStack(itr,:) = absVel(sIdx:eIdx);
            dspStack(itr,:) = absDsp(sIdx:eIdx);
           
            if ( ismember(itr,hundreds) & opts.verbose )
                gcf; plot(idx_stack,log10(accStack(itr,:)),'k','lineWidth',2);
                ylm = get(gca,'ylim');
                line([0 0],ylm,'color','r','lineWidth',2)
                pause(.5)
            end
            
        else
            fprintf(1,'Skipping trace, not enough samples around pick\n')
            flg_shortSignal(itr) = true;
        end
    end
end

% Package output
stacks.acc = sum(accStack);
stacks.vel = sum(velStack);
stacks.dsp = sum(dspStack);

counts.wrongSR       = sum(flg_wrongSR);
counts.wformNotFound = sum(flg_wformNotFound);
counts.shortSignal   = sum(flg_shortSignal);
