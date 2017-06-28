function [ctlg] = import_NCAeqDD(fileFullName)

%% Import data from text file.
% Script for importing data from the following text file:
%
%    /scratch/memeier/data/nocal/waldhauser/NCAeqDD.v201112.1.txt
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2016/02/02 20:07:28

% For more information, see the TEXTSCAN documentation.
formatSpec = '%4f%3f%3f%3f%3f%7f%12f%13f%9f%8f%8f%4f%8f%5f%f%[^\n\r]';
fileID     = fopen(fileFullName,'r');

% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);

% Create output variable
ctlg = table(dataArray{1:end-1}, 'VariableNames', {'yr','mt','dy','hr','min','sec','lat','lon','dep','eh1','eh2','az','ez','m','id'});

clearvars fileFullName formatSpec fileID dataArray ans;