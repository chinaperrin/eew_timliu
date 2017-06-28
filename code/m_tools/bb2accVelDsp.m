function [acc,vel,dsp,dsp2] = bb2accVelDsp(vh,sr,ppxIdx,intMode,fHighpass,fOrder,fMode)

o.intMode = intMode;

vhd   = diff(vh)*sr;                                                       % Differentiate to acc (from [m/s] to [m/s^2])
vhd   = [vhd; vhd(end)];                                                   % Duplicate last sample to have equal length
vhi   = integrate_wform(vh,sr,ppxIdx,o);                             % Integrate to dsp (from [m/s] to [m])
vhih  = bworth(vhi,sr,fHighpass,'high',fOrder,fMode);                      % Highpass

acc  = vhd;
vel  = vh;
dsp  = vhih;
dsp2 = vhi;   % Highpassing dsp trace removes all static offset