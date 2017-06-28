function [out] = read_sac_trace3_bigEndian(traceFullName)

% Read sac-trace with RSAC.m
%
% If waveform is not found, returns empty arrays
% If waveform is not ok, returns the read values, but raises flgIssue
% menandrin@gmail.com, 130304

global o
flgSkip = false;

% Get non-available waveforms from bigstar?
if ~exist(traceFullName,'file')
    if o.scp_wforms
        scp_wform(traceFullName)
    else
        fprintf(1,'Waveform not found. Switch o.scp_wforms on if you want to download it from bigstar\n')
        flgSkip = true;
    end
end

if ~flgSkip
    % Read data
    S    = rsac_bigEndian(traceFullName);
    sraw = S.amplitudes;
    t    = S.timevector';
    dt   = roundn(S.header(1),-4);
    sr   = 1/dt;
    
    if (sr>999)
        fprintf(1,'8UNG: rounding is not appropriate for such sampling rates.\n')
    end
    
    % Absolute reference time
    year = S.header(71);
    jday = S.header(72);
    hr   = S.header(73);
    mm   = S.header(74);
    sec  = S.header(75);
    msec = S.header(76);
    
    out.stla = S.header(32);
    out.stlo = S.header(33);
    out.stel = S.header(34);
    out.stdp = S.header(35);
    
    mdate_ref   = datenum(year-1, 12, 31, hr, mm, sec) + jday + msec/86400/1e3;     % Reference time in [days]
    tref        = datestr(mdate_ref,'yyyy/mm/dd HH:MM:SS.FFF');                     % Absolute reference time
    
    dt          = S.header(6);                                                      % dt between ref time and 
    mdate_start = mdate_ref+ dt/86400;                                              % start of recording [sec]
    tstart      = datestr(mdate_start,'yyyy/mm/dd HH:MM:SS.FFF');                   % Absolute reference time
    
    
    % Test for issues
    ns = numel(sraw);
    nt = numel(t);
    if ~isequal(ns,nt)
        flgIssue = true;
        fprintf(1,'\n\t\tns ~= nt, something wrong with waveform, skip.\n')
    else
        flgIssue = false;
    end
    
    % Write to out structure
    out.sraw     = sraw;
    out.t        = t;
    out.sr       = sr;
    out.tstart   = tstart;
    out.tref     = tref;
    out.flgIssue = flgIssue;
else
    out.flgIssue = true;
end