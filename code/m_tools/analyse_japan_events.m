clear all

addpath(genpath('../m_plot/'))
addpath(genpath('../../../../../matlab/fct_downloads/'))
addpath(genpath('../../../../../matlab/myFct/src/'))


% NEIC event list for magnitudes
neicListFullName = '/scratch/memeier/VS/data/wform/japan/neic/m5p_since1996_zmax35.txt';
neicList         = import_neicList(neicListFullName);

% Get a list of sub-directories (=events) in the specified wform-directory
wformDir = '/scratch/memeier/VS/data/wform/japan/k_kik/M6p_D25km/';
dirList  = dir(wformDir);
eqIdx    = cellfun(@(x) str_isnumeric(x), {dirList.name}');
eqList   = dirList(logical(eqIdx));
neq      = size(eqList,1);
fprintf(1,[' ',num2str(neq),' waveform directories / events:\n\n'])

dt      = zeros(neq,1);
idxNeic = zeros(neq,1);

% For each sub-directory (= event) ...
for ieq = 1:neq
    
    eqName     = getfield(eqList,{ieq}, 'name');
    eqFullName = strcat(wformDir,eqName);
    fprintf(1,[' ... event directory ',num2str(ieq),' out of ',num2str(neq),': ',eqName,'\n'])
    
    % Find all k-Net waveform files with correct extension
    fileList_knt = dir(strcat([eqFullName,'/knt/']));   % get a list of all files
    fileList_knt = {fileList_knt.name}';
    wfIdx_knt    = regexp(fileList_knt,'[ENU][WSD]$');
    isWf_knt     = cellfun(@(x) ~isempty(x), wfIdx_knt);
    wfList_knt   = fileList_knt(isWf_knt);
    nwf_knt      = numel(wfList_knt);
    someFileName = [eqFullName,'/knt/',wfList_knt{1}];
    
    fileNameTxt = strcat(['/scratch/memeier/VS/data/wform/japan/neic/maps/',eqName,'.txt']);
    diary(fileNameTxt)
    [neic_idx,differences] = find_neic_list_entry(eqName,someFileName,neicList);
    diary off
    
    % Print map
    fileName = strcat(['/scratch/memeier/VS/data/wform/japan/neic/maps/',eqName,'.png']);
    gcf;
    set(gcf,'PaperPositionMode','auto') % --> ensures that sizes are the same as on screen
    print('-dpng',fileName)
    
    1+1;
end