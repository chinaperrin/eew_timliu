function [indices] = find_station_traces(stationName,trList)

matches = regexp(trList.fullName,stationName,'match');
indices = find(cellfun(@(x) ~isempty(x),matches));