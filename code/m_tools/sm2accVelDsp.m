function [acc,vel,dsp,dsp2] = sm2accVelDsp(ah,sr,ppxIdx,intMode,fHighpass,fOrder,fMode)

o.intMode = intMode; 

% ah = highpass filtered input accelerogram
ahi     = integrate_wform(ah   ,sr,ppxIdx,o);                        % Integrate to vel (from [m/s^2] to [m/s])
ahih    = bworth         (ahi  ,sr,fHighpass,'high',fOrder,fMode);   % Highpass
ahihi   = integrate_wform(ahih ,sr,ppxIdx,o);                        % Integrate to dsp (from [m/s] to [m])
ahihih  = bworth         (ahihi,sr,fHighpass,'high',fOrder,fMode);   % Highpass

acc  = ah;                                                  
vel  = ahih;
dsp  = ahihih;
dsp2 = ahihi;   % Highpassing dsp trace removes all static offset