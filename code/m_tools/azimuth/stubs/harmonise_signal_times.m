function [Z,N,E,t,pxIdx] = harmonise_signal_times(z,n,e,tz,te,tn,pxIdx)

% Harmonise vector lengths
commonStartTime = max([tz(1),tn(1),te(1)]);
commonEndTime   = min([tz(end),tn(end),te(end)]);

% Vertical
[~,startIdx_Z] = min(abs(tz-commonStartTime));
[~,endIdx_Z]   = min(abs(tz-commonEndTime));
t              = tz(startIdx_Z:endIdx_Z);
Z              = z (startIdx_Z:endIdx_Z);
pxIdx          = pxIdx-startIdx_Z;

% North
[~,startIdx_N] = min(abs(tn-commonStartTime));
[~,endIdx_N]   = min(abs(tn-commonEndTime));
N              = n(startIdx_N:endIdx_N);

% East
[~,startIdx_E] = min(abs(te-commonStartTime));
[~,endIdx_E]   = min(abs(te-commonEndTime));
E              = e(startIdx_E:endIdx_E);