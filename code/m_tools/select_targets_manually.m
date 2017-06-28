% Criteria for target

% M6.5 for EGU presentation
%idxList = find((zList.m>6.3) & (zList.m<7) & (zList.hypDist<15) & (zList.hypDist>10));
%idx = idxList(7);       % Select one
% --> these are indices wrt/ traceLists, e.g. zList

% M3.0 for EGU presentation
idxList = find((zList.m>2.9) & (zList.m<3.2) & (zList.hypDist<35) & (zList.hypDist>30));
idx = idxList(217);       % Select one

eventId = eqs.eventIdTrList(idx);       % eventId
eqsIdx  = find(eqs.eventId==eventId);   % eqs-index of this event Id
% eqsIdx = find(eqs.eventId==eqs.eventIdTrList(idx)); would be the more
% direct form for finding eqsIdx

% Load regression coefficients
C = Ca2{eqsIdx};

% Training set
idxTarget = eqs.traceId{eqsIdx};
ntarget   = numel(idxTarget);

% Check trace names
zList.fullName(idx)
zList.fullName(idxTarget)

idxTraining = setdiff(1:nz,idxTarget); % traceList index of all other traces
zTraining   = zList.selectSubList(idxTraining);
eTraining   = eList.selectSubList(idxTraining);
nTraining   = nList.selectSubList(idxTraining);
wt          = Wt(idxTraining);

% Target set
zTarget = zList.selectSubList(idx);
eTarget = eList.selectSubList(idx);
nTarget = nList.selectSubList(idx);

zFullName = zTarget.fullName{1};
eFullName = eTarget.fullName{1};
nFullName = nTarget.fullName{1};
fprintf(1,'\t%s\n','Selected traces:',zFullName,eFullName,nFullName)
fprintf(1,'\t%s\n','Other traces from same event:',zList.fullName{idxTarget})

if ( (~exist(zFullName,'file')) && (o.scp_wforms) )
    slash_idx        = regexp(zFullName,'/');
    zTraceName       = zFullName(slash_idx(end)+1:end);
    pathName         = zFullName(1:slash_idx(end));
    [recordName]     = get_recordName(zFullName);
    threeCompPattern = strcat([pathName,recordName,'*']);
    
    scp_wform(threeCompPattern)
    %cp_wform_lacie(threeCompPattern,pathName)
end

% Gather meta-info
stName = zList.stationName{idx};
stLat  = zList.stationLat(idx);
stLon  = zList.stationLon(idx);
eqLat  = zList.eqLat(idx);
eqLon  = zList.eqLon(idx);
eqZ    = zList.eqZ(idx);
eqDate = zList.eqDate{idx};

% Save catalog values
target.m_ctlg = zTarget.m;
target.r_ctlg = zTarget.hypDist;
target.stName = stName;