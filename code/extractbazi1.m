%This program extracts the back azimuth estimates and true back azimuth
%from the TraceList. The function must be in the same directory as
%the script "load_traceList_data." The function outputs text files with
%a list of the back azimuth estimates. The function extracts from a 
%single back azimuth estimate technique, which must be changed in the
%code. This code is derived from extractbazi.m
%
%Author:        Timothy Liu
%Last Revised:  6/21/17

%clear everything out
clear all
%load the traceList
load_traceList_data;

%names of techniques for estimating back azimuth
%ename = 'trueval';
%ename = 'onsite';
%ename = 'sac';
%ename = 'deichmann';
ename = 'la05';

%set LEN to 15 for any technique except trueval, which should be 1
LEN = 15;

%iterate through estimates at each time step
for i = 1:LEN
    %name file with technique name and time step
    filename = char(strcat(ename, '_', int2str(i), '.txt'));
    fileID = fopen(filename, 'w');
    fprintf('Outputting %s\n', filename);
    %print records for all event_triggers - nz is number of triggered
    %stations summed across all events
    for x = 1:nz
        %uncomment whichever estimate technique is being used
        
        %baze = zList.scalFeature{x}.bazi.trueVal(i);
        %baze = zList.scalFeature{x}.bazi.onsite.baziHat(i);
        %baze = zList.scalFeature{x}.bazi.sac.baziHat(i);
        %baze = zList.scalFeature{x}.bazi.deichmann.baziHat(i);
        baze = zList.scalFeature.bazi{x}.la05.baziHat(i);
        fprintf(fileID, '%.4f\n', baze);
    end
    fclose(fileID);
end
