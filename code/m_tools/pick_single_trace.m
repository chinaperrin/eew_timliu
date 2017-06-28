clear all

traceFullName = '/scratch/memeier/data/socal/scsn_040101_051231/cms/10147909/10147909.AZ.SOL.HHZ.sac';

ntap = 100;
npad = 100;

% Picker & SNR
prePxWindow = 5;              % No. of sec before theoretical pick that pick search is started
snr_min     = 20;             % Minimum SNR as returned by the picker from the vertical records
fLow_px     = 3;
fUp_px      = 10;
pxThreshold = 5;
sWindow     = 1;              % Window lenght for measuring signal amps after p-pick [sec]
nWindow     = 1;              % Window lenght for measuring noise before p-pick [sec]
wmax        = 1e1;
wmin        = .1;
sWtFct      = 'lin';
nWtFct      = 'lin';



[sraw,t,sr,flgIssue] = read_sac_trace2(traceFullName);
s                    = sraw/100;    % converting from [cm/s]   to [m/s]
dt                   = 1/sr;    
ts                   = t(1);


smpx      = sraw - mean(sraw(1:2*ntap));
stappx    = taper(smpx,sr,ntap);
[spx,~,~] = butter_pass_tdomain_f(stappx,fLow_px,fUp_px,sr,8,npad,0);


[ppxIdx_z,snr_z,~,~] = SBPx(spx,1/sr,0,pxThreshold,sWindow,nWindow,wmax,wmin,sWtFct,nWtFct,ntap,1);
