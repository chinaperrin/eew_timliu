function out = upscale_presignal_noise(s,ppxIdx,sr,p95target,opts)
% Takes a signal, cuts out part of the pre-signal noise, replicates and
% concatenates the noise, upscales it, then adds the upscaled noise to the
% original signal; this summation is repeted nrnd times, with random onset
% points. 
% 
% The function is a take-out from the more comprehensive function
% recompute_Pd_with_added_noise.m

% Unpack opts
nrnd      = opts.nrnd;

% Other parameters
nzcMin = 10;
ns     = numel(s);
ntap   = 100;

%opts.debug=0;

out.wforms      = cell (nrnd,1);
out.scaleFactor = zeros(nrnd,1);

% Cut out 'reliable' noise part
if     ppxIdx>2000; isn = 1000;          % If there is a lot of pre-signal waveform, cut out more...
elseif ppxIdx>1000; isn = 500;
else                isn = 2*ntap;
end
ien   = round(.95*ppxIdx);
snRaw = s(isn:ien);
snRaw = snRaw - mean(snRaw);

if opts.debug; clf; hold on; plot(s,'y'); plot(isn:ien,snRaw,'r'); set(gca,'xlim',[1,ppxIdx+40]); end

% Concatenate replications of snRaw: i) find last three zero crossings,
% ii) find last positive maximum in between these, iii) do same at
% beginning of snRaw, iv) use the segment between those indices
% (=snSegment) for concatenation.
if opts.debug; oc.plot = true; else oc.plot=false; end
crIdx = get_zero_crossings(snRaw,oc);
nzc   = numel(crIdx);

if nzc>=nzcMin
    [~,iFirstMax]   = max(snRaw(1:crIdx(2)+1));
    [~,iLastMaxRel] = max(snRaw(crIdx(end-2):crIdx(end)));
    iLastMax        = crIdx(end-2)+iLastMaxRel-1;
    snSegment       = snRaw(iFirstMax:iLastMax);
    nsseg           = numel(snSegment);
    
    % Replicate snSegment to be long enough for random onset index sampling
    nspx        = ppxIdx+5*sr;   % Length of signal to which noise is added
    if nspx>ns; nspx = ns; end
    nscatTarget = 5*nspx;        % Make length of noise signal much longer
    nrep        = ceil(nscatTarget/nsseg);
    snCat       = [snRaw(1:iLastMax); repmat(snSegment,nrep,1)];
    idxCat      = [iLastMax;          iLastMax+(1:nrep)'*nsseg];
    nscat       = numel(snCat);
    if opts.debug; clf; hold on; plot(snCat,'y'); plot(snRaw,'r'); plot(idxCat,0,'dm','markerSize',12); end
    
    % Scale concatenated noise vector snCat to target amplitude
    p95      = prctile(snCat,95);
    scFactor = p95target/p95;
    if scFactor<1; scFactor = 0; end    % if noise level is higher than target, don't add any noise.
    snScaled = scFactor*snCat;
    if opts.debug; gcf; hold on; plot(snScaled,'w'); end
    
    % Choose random starting point and add upscaled noise to signal
    for irnd = 1:nrnd
        
        %print_iteration_numbers(irnd,nrnd,'hundreds')
        
        imax   = nscat-nspx-1;
        if imax<1; imax=1; end
        isrand = randi(imax,1,1);
        ierand = isrand+nspx-1;
        if ierand>nscat; ierand=nscat; end
        noise    = snScaled(isrand:ierand);
        sOrg     = s(1:nspx);
        sNoised  = sOrg+noise;
        p95check = prctile(sNoised(ntap:round(.95*ppxIdx)),95);         % Check p95 of noise after summation
        
        out.wforms{irnd}      = sNoised;
        out.scaleFactor(irnd) = scFactor;

        if opts.debug;
            clf; hold on; plot(sOrg,'r'); plot(noise,'y'); plot(sNoised,'m')
            fprintf(1,sprintf('p95: %em/s/s before, %em/s/s after scaling, target=%em/s/s\n',p95,p95check,p95target))
        end
    end
else
    fprintf(1,'Not enough zero crossings...\n')
end