clear all
addpath(genpath('m_tools/'))

% Import data set batches
dataSetNames = {'../data/socal/M4/out/i39/causal_0p5_gba/'; ...
                '../data/socal/M5/out/i39/causal_0p5_gba/'; ...
                '../data/socal/M6p/out/i39/causal_0p5_gba/'};

TraceList = import_trLists(dataSetNames,'traceList');


% Separate entries for vertical and horizontal records
fprintf(1,'\nSplitting TraceList into z-, e- & nList ... ')
idx_Z = find(strcmp(TraceList.station.ocode,'Z'));
idx_E = find(strcmp(TraceList.station.ocode,'E'));
idx_N = find(strcmp(TraceList.station.ocode,'N'));

zList = TraceList.selectSubList(idx_Z);
eList = TraceList.selectSubList(idx_E);
nList = TraceList.selectSubList(idx_N);
%clear TraceList
 
nz = numel(zList.eq.m);
ne = numel(eList.eq.m);
nn = numel(nList.eq.m);


% Associate records with events
eqs = compile_eqs_structure(zList,[]);

% Time interval at which back-azimuth ("bazi") is estimated; in sec since
% P-wave onset
tbazi = zList.prop.azi.bazIntervalVector;

% For the 45th record, the estimated bazi from method "onsite" is stored as
% follows:
zList.scalFeature{45}.bazi.onsite.baziHat
% >> a 15-by-1 vector, containing the estimates using the 15
%    different signal lengths stored in tbazi

% Getting the estimates for all the entries is a little complicated (sorry): 
baziHat = cell2mat(cellfun(@(x) x.bazi.onsite.baziHat', zList.scalFeature,'uniformOutput',0));

% The true back-azimuth is stored here:
trueBazi = zList.scalFeature{45}.bazi.trueVal;
baziTrue = cell2mat(cellfun(@(x) x.bazi.trueVal', zList.scalFeature,'uniformOutput',0));
