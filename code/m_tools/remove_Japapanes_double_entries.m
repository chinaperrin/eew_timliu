clear all

doubleDirList = importdata('/scratch/memeier/data/japan/k_kik/M5p_z25km/doubleDirs.txt');
ndb           = numel(doubleDirList);

% load M5 & M4 trLists
clear m4List
m5ListFullName = '/scratch/memeier/data/japan/k_kik/M5p_z25km/out/i38/causal_0p5/trList.mat';
%m4ListFullName = '/scratch/memeier/data/japan/k_kik/M4_z25km_96to03/out/i38/causal_0p5/trList.mat';
%m4ListFullName = '/scratch/memeier/data/japan/k_kik/M4_z25km_04to10/out/i38/causal_0p5/trList.mat';
m4ListFullName = '/scratch/memeier/data/japan/k_kik/M4_z25km_11to15/out/i38/causal_0p5/trList.mat';
load(m5ListFullName)
m5List = TraceList;
load(m4ListFullName)
m4List = TraceList;
n4     = numel(m4List.m);

% for each event in doubleDirList, find entries in both lists, check if
% they are really the same, remove them from M4-list and re-save M4-list
keepMe = true(n4,1);
for idb = 1:ndb
    dirLongName = doubleDirList{idb};
    slashIdx    = regexp(dirLongName,'/');
	dirName     = dirLongName(slashIdx(end)+1:end);
      
    idx4 = find(cellfun(@(x) ~isempty(x), regexp(m4List.fullName,dirName)));
    %idx5 = find(cellfun(@(x) ~isempty(x), regexp(m5List.fullName,dirName)))
    keepMe(idx4)=false; 
end
sum(keepMe==0)
idxKeep = find(keepMe);
m4newList = m4List.selectSubList(idxKeep);
TraceList=m4newList;
%save(m4ListFullName,'m4newList')


%m4ListFullName = '/scratch/memeier/data/japan/k_kik/M4_z25km_96to03/out/i38/causal_0p5/trList.mat';
%m4ListFullName = '/scratch/memeier/data/japan/k_kik/M4_z25km_04to10/out/i38/causal_0p5/trList.mat';
%m4ListFullName = '/scratch/memeier/data/japan/k_kik/M4_z25km_11to15/out/i38/causal_0p5/trList.mat';
%load(m4ListFullName)
%TraceList=m4newList;
%save(m4ListFullName,'TraceList')
%clear TraceList m4newList


% %% Remove Japanes doubly contained traces
% ntr = numel(TraceList.m);
% idxDoubleList = [];
% for itr =8001:ntr
%     dsn    = TraceList.dataSetName{itr};
%     isJapa = strcmp(dsn,'kNet') |strcmp(dsn,'kikNet');
%     m      = TraceList.m(itr);
%     if m>=5 &isJapa
%         %fullName = 'M4_z25km_04to10/20060909024000/kik/CHBH040609090240.UD2';
%         fullName  = TraceList.fullName{itr};
%         slashIdx  = regexp(fullName,'/');
%         traceName = fullName(slashIdx(end)+1:end);
%         isInM5Dir = ~isempty(regexp(fullName,'M5p_z25km'));
%         
%         % Check if it also has an entry in one of the M4 directories
%         if isInM5Dir
%             
%             idx = find(cellfun(@(x) ~isempty(x), regexp(TraceList.fullName,traceName)));
%             
%             if numel(idx)>1
%                 % Find out which one(s) are from an M4-dir; block them
%                 for ii = 1:numel(idx)
%                     thisFullname = TraceList.fullName{idx(ii)};
%                     if ~isempty(regexp(thisFullname,'M4_z25km'))
%                         idxDoubleList = [idxDoubleList; idx(ii)];
%                     end
%                 end
%             end
%         end
% %         searchPattern = sprintf('*M4/',) 
% %         regexp(fullName,'M4*')
%         
%     end
% end
