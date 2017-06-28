%This program extracts the back azimuth estimates and true back azimuth
%from the TraceList. The function must be in the same directory as
%the script "load_traceList_data." The function outputs text files with
%a list of the back azimuth estimates
%
%Author:        Timothy Liu
%Last Revised:  6/21/17

%clear everything out
clear all
%load the traceList
load_traceList_data;
x = 1;

%names of techniques for estimating back azimuth
ename = {'trueval', 'onsite', 'sac', 'deichmann', 'la05'};
%array of back azimuth estimates - each element except the first
%is a list of 15 estimates made at 15 different times steps
%first element is the true azimuth value
bazes = {zList.scalFeature{x}.bazi.trueVal...
    zList.scalFeature{x}.bazi.onsite.baziHat, ...
    zList.scalFeature{x}.bazi.sac.baziHat,...
    zList.scalFeature{x}.bazi.deichmann.baziHat,...
    zList.scalFeature{x}.bazi.la05.baziHat};
%iterate through back azimuth estimate techniques and print to txt file
for tech = 1:length(bazes)
    %object with estimates from one technique at different times
    baz_t = bazes(tech);
    %iterate through estimates at each time step
    for i = 1:length(baz_t{1,1})
        %name file with technique name and time step
        filename = char(strcat(ename(tech), int2str(i), '.txt'));
        fileID = fopen(filename, 'w');
        fprintf('Outputting %s\n', filename);
        %print records for all event_triggers - nz is number of triggered
        %stations summed across all events
        for x = 1:nz
            fprintf(fileID, '%.4f\n', baz_t{1, 1}(i));
        end
        fclose(fileID);
    end
end