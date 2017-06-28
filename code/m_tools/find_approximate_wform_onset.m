function [idx] = find_approximate_wform_onset(s,dt,pctThreshold)

o.plot = false;
%pctThreshold  = 0.001;
%traceFullName = '/scratch/memeier/data/ngawest/wform/CHICHI00_TCU072_999.AT2';
%[sraw,meta]   = read_ascii_wform_nga(traceFullName,1);
%ns            = numel(sraw);
%s             = sraw*9.80665;   % From [g] to [ms^-2] (BIPM: 1g = 9.80665 ms^-2)
%sr            = meta.sr;    
%dt            = 1/sr;

cav     = cumsum(abs(s))*dt;
cavn    = cav/cav(end);
[~,idx] = min(abs(cavn-pctThreshold));

if o.plot
    clf; hold on; box on; grid on;
    plot(cavn)
    plot(s,'r')
    plot(idx,0,'xr')
end

%cavnh = bworth(cavn,sr,.1,'high',fOrder,'causal');   % Bandpass filter for picking