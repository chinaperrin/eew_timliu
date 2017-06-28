function idxThd = get_threshold_indices(tseries,thresholdList)

nth    = numel(thresholdList);
idxThd = zeros(nth,1);
for ith = 1:nth
    idx = find(abs(tseries)>=thresholdList(ith),1,'first');
    if ~isempty(idx); idxThd(ith)=idx; end
end
%idxThd(cellfun(@isempty, idxThd)) = {0};
%idxThd                            = cell2mat(idxThd);