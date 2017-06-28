% HOW TO USE TRACELIST OBJECTS
% 
% TraceLists are a way to organise waveform meta-data in matlab. They are
% matlab objects that contain a list-entry for each record in a waveform
% data set. Each entry has various attributes such as station meta-
% information, the waveform file name and its path, computed waveform 
% features, etc. They basically contain all the information except the 
% waveform itself. The main advantages of traceLists are i) that data from 
% various repositories are merged into a single ready-to-use data format, 
% and ii) the simplicity with which data sub-sets can then be handled. 
% This script shows how traceLists can be used and manipulated.
%
% Note: to load a traceLists you need to be in the directory where the 
% object definition file ("traceList.m") is located, or add the corresonding
% path, e.g. with "addpath(genpath('../../../path-to-traceList.m/'))"
%
% All ground motion values are in SI units.
%
% menandrin@gmail.com, 170620

clear all

% Load 1 traceList each for vertical and two horizontal components; note
% that the two horizontal components are not always N-S and E-W 
load('../data/socal/M5/out/i39/causal_0p5_gba/trList.mat')

% Number of records in list:
ntr = numel(TraceList.eq.m);

% There are a lot of attributes for each trace (not all of which are populated) ... 
TraceList

TraceList.station.name    % Name of recording station
TraceList.station.lat     % Station latitutde 
TraceList.station.ocode   % Orientation of sensor component
TraceList.eq.m            % Earthquake magnitude
TraceList.eq.lat          % Earthquake latitude
TraceList.eq.date         % Earthquake origin time
TraceList.eq.t0           % 
TraceList.dist.hyp        % Hypocentral distance
TraceList.fullName        % Name of waveform file including path
TraceList.px.p.t          % P-wave onset time, in seconds since origin time
TraceList.pga.pns.ampi    % Peak observed acceleration in [m/s/s]
TraceList.var.v1          % v1-v8 are fields of cell type, to store pretty 
                          % much anything that is not covered by the existing fields

% So you can plot, for example, 
clf; plot(TraceList.dist.hyp,TraceList.px.p.t,'x');

% In summary: there is a ton of information in these traceLists. If you
% need something, just ask me. Chances are I can just point you to a field
% in the traceLists.





%% Select a subList
%  The "selectSubList" function allows to extract a sub-set of the data,
%  e.g. all records from events with a certain magnitude. 
idx      = find(TraceList.eq.m>5.7);
testList = TraceList.selectSubList(idx);


% You can use it to make a separate traceLists for each sensor component 
% (vertical, east, north) ...
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
fprintf(1,' done.\n')

% ... or to extract records with hypocentral distances <25km
idx    = find(zList.dist.hyp<25);
nsList = zList.selectSubList(idx);

% ... extract records with hypocentral distances <25km & magnitudes>5.4  
tmpList = zList.selectSubList(find(zList.dist.hyp<50 &zList.eq.m>5.4));


% --> this selection function is one of the most useful features of the 
%     traceList format. You make a selection of records and you 
%     automatically select all attributes for those records. I think that 
%     was why I started it in the first place :-)





%% Various commands
%  Print main info on 51th trace
itr=51; % --> traceList index of randomly chosen record
zList.printSingleTraceSummary(itr);

% The three lists have the same order, i.e. you can use the same index to
% find the three components of the same record in the three lists:
zList.printSingleTraceSummary(itr);
eList.printSingleTraceSummary(itr);
nList.printSingleTraceSummary(itr);

% The size of the object
zList.printObjectSize;

% zList.prop contains a couple of parameters that were used for generating 
% the traceList ...
zList.prop

% ... e.g. the frequency bands that are used by the Gutenberg Algorithm
zList.prop.fc






%% Add two lists together
% e.g. make list of all records of the two large Turkish events concatenate them.
idx1            = find(strcmp(zList.eq.name,'3147406'));
firstEventList  = zList.selectSubList(idx1);
idx2            = find(strcmp(zList.eq.name,'7062511'));
secondEventList = zList.selectSubList(idx2);

bothEventsList = traceList(0);                  % Intitiate traceList
bothEventsList.appendList(firstEventList);      % Append first list
bothEventsList.appendList(secondEventList);     % Append second list
% NOT SURE WHY THIS DOESNT WORK. NEED FIX.



%% Clone list
%  Because traceLists are "handle classes" copying them with a command like
%  newzList = zList only creates a new reference to the same object, not a
%  new object itself. If you then make changes to newzList the same changes
%  will also be applied to zList itself. To create a new independent
%  instance of a traceList, "clone" the list using the selectSubList function:
idxAll              = 1:numel(zList.eq.m);
newIndependentzList = zList.selectSubList(idxAll);



%% Remove unneeded fields to save space
%  This can be done using the following external function, specifying the 
%  fields that you would like to discard
overwrite_unneeded_fields(zList,{'var.v7';'fb.cav';'fb.amax'})






% If you've made it that far, olé! I'm always happy for comments and
% suggestions.
%
%
%                             I liked ascii-tables better...
%             ,;::\::\        
%           ,'/' `/'`/       /
%       _\,: '.,-'.-':.     /
%      -./"'  :    :  :\/, /
%       ::.  ,:____;__; :-       
%       :"  ( .`-*'o*',);       
%        \.. ` `---'`' /
%         `:._..-   _.'
%         ,;  .     `.
%        /"'| |       \
%       ::. ) :        :
%       |" (   \       |
%       :.(_,  :       ;
%        \'`-'_/      /
%         `...   , _,'
%          |,|  : |
%          |`|  | |
%          |,|  | |
%      ,--.;`|  | '..--.
%     /;' "' ;  '..--. ))
%     \:.___(___   ) ))'