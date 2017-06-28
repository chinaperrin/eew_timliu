function [countList] = count_station_instances(staList)
% Go through list of station and count instances

n         = numel(staList);
countList = cell(n,2);
c_sta     = 0;

while (~isempty(staList))
   
    sta_ct = staList{1};
   
    % If station has not been read before, add it to countList and count
    % instances, else skip it 
   
    isitahit = regexp(sta_ct,countList(:,1));
    idx      = find(cellfun(@(x) ~isempty(x), isitahit));

    if (isempty(idx))
   
        c_sta      = c_sta + 1;
        
        isitahit   = regexp(staList,sta_ct);
        idx        = find(cellfun(@(x) ~isempty(x), isitahit));
        ninstances = numel(idx);
        
        countList{c_sta,1} = sta_ct;
        countList{c_sta,2} = ninstances;
        
        % Throw out found list entries (to speed things up)
        lgc      = true(numel(staList),1);
        lgc(idx) = false;
        staList  = staList(lgc);
        
    end
end

% Clear out unfilled entries
last_idx  = find(cellfun(@(x) ~isempty(x), countList(:,1)),1,'last');
countList = countList(1:last_idx,:);

% Sort list
[~,idx_std] = sort(cell2mat(countList(:,2)),'descend');

countList(:,1) = countList(idx_std,1);
countList(:,2) = countList(idx_std,2);

