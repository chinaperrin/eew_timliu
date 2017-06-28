function [accNoiseAmps,velNoiseAmps,dspNoiseAmps] = get_absolute_noise_amps(trList,mRange)

idx       = find(logical(trList.m>=mRange(1)) & logical(trList.m<mRange(2)));
shortList = trList.selectSubList(idx);

n            = numel(shortList.m);
nvalsMax     = max(cellfun(@(x) numel(x), shortList.accNoise));
accNoiseAmps = zeros(n,nvalsMax);
velNoiseAmps = zeros(n,nvalsMax);
dspNoiseAmps = zeros(n,nvalsMax);

for i = 1:n

    acc   = cell2mat(shortList.accNoise(i))';
    vel   = cell2mat(shortList.velNoise(i))';
    dsp   = cell2mat(shortList.dspNoise(i))';
    nvals = numel(acc);
    
    % If waveform wasn't long enough, fewer noise amplitude samples are
    % available. Pad with zeros.
    if nvals~=nvalsMax
        acc = [zeros(1,nvalsMax-nvals), acc];
        vel = [zeros(1,nvalsMax-nvals), vel];
        dsp = [zeros(1,nvalsMax-nvals), dsp];
    end
    accNoiseAmps(i,:) = acc;
    velNoiseAmps(i,:) = vel;
    dspNoiseAmps(i,:) = dsp;
end
