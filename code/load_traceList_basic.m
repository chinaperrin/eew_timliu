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