function [indices] = find_event_traces(eventName,trList)

matches = regexp(trList.fullName,eventName,'match');
indices = find(cellfun(@(x) ~isempty(x),matches));