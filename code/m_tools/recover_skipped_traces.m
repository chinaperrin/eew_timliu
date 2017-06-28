clear all

global iN

iN = 33;
fMode           = 'causal';      % 'causal'/'acausal'
outDirAppendix = '_0p5';
outDirName     = sprintf('%s%s',fMode,outDirAppendix);

%dataSetBaseNames = {'/scratch/memeier/data/socal/scsn_900101_011231/out/'};
%nds              = numel(dataSetBaseNames);
        
dataSetBaseNames = {'/scratch/memeier/data/ngawest/out/'; ...
                    '/scratch/memeier/data/japan/k_kik/M5p_z25km_surface/out/'; ...
                    '/scratch/memeier/data/socal/scsn_900101_011231/out/'; ...
                    '/scratch/memeier/data/socal/scsn_020101_031231/out/'; ...
                    '/scratch/memeier/data/socal/scsn_040101_051231/out/'; ...
                    '/scratch/memeier/data/socal/scsn_060101_071231/out/'; ...
                    '/scratch/memeier/data/socal/scsn_080101_091231/out/'; ...
                    '/scratch/memeier/data/socal/scsn_100101_111231/out/'; ...
                    '/scratch/memeier/data/socal/scsn_120101_131231/out/'; ...
                    '/scratch/memeier/data/socal/scsn_140101_151231/out/'; ...
                    '/scratch/memeier/data/socal/scsn_120101_151231_M2_2p5/out/'; ...
                    '/scratch/memeier/data/socal/scsn_120101_151231_M2p5_3/out/'};

                
dataSetNames             = strrep(dataSetBaseNames,'/out/',sprintf('/out/i%i/%s/',iN,outDirName));                
[SkipList ,fc,snpLength] = import_trLists(dataSetNames,'skipList');
[TraceList,fc,snpLength] = import_trLists(dataSetNames,'traceList');

summary_trList(TraceList,SkipList,1,0)

rmax = 30;
mmin = 6;

idx_NS = find(SkipList.hypDist<=rmax &SkipList.m>=mmin);
nsList = SkipList.selectSubList(idx_NS);
unique(nsList.comment)
summary_trList([],nsList,1,0)


fprintf(1,'\nSplitting TraceList into z-, e- & nList ... ')
idx_Z = find(strcmp(nsList.orntCode,'Z'));
idx_E = find(strcmp(nsList.orntCode,'E'));
idx_N = find(strcmp(nsList.orntCode,'N'));
znsList = nsList.selectSubList(idx_Z);
ensList = nsList.selectSubList(idx_E);
nnsList = nsList.selectSubList(idx_N);
nz    = numel(znsList.m);
ne    = numel(ensList.m);
nn    = numel(nnsList.m);

unique(nsList.comment)
max(znsList.m)

ntr = numel(znsList.m);
for itr = 1:ntr

end
%    'Z and H traces have very different lengths'
%     'clipped'
%     'has unclipped BB record'
%     'less than 3 corecs'
%     'low sr'
%     'low vertical SNR'
%     'more than one sr'
%     'most time-values negative'
%     'no StaLtaPx'
%     'not enough data for SNR'
%     'outlier'
%     'too far'
%     'too many identical amps'
%     'waveform too long'
%     'waveform too short'
%     'wrong pick'