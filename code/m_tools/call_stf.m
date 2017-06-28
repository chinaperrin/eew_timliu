function [stf,t,stdout] = call_stf(M0,sd,pv,rv,th,k1,k2,sr,ns,ne)
% This is a warpper to call stf.f from Nicholas Deichmann.
%
% INPUT:    M0
%           sd
%           pv
%           rv
%           th
%           k1
%           k2
%           sr
%           ns
%           ne

%addpath(genpath('~/programs/colleagues/nicolasDeichmann/stf/'))
exePath         = '~/programs/colleagues/nicolasDeichmann/stf/';
outFileFullName = 'stf.dat';
%outFileFullName = sprintf('%sstf.dat',exePath);
if exist(outFileFullName,'file'); unix(['rm ', outFileFullName]); end


% 1. Convert M0 to exponent & mantissa
m0exponent = floor(log10(M0));
m0mantissa = M0/10^m0exponent;
%test = m0mantissa*10^m0exponent; isequal(M0,test)

% Floating point precision for mantissa = D. This is necessesary because
% stf.f requires integer input
D = 8;
m0mantissaD = round(10^8*m0mantissa);
m0exponentD = m0exponent-D;
%test = m0mantissaD*10^m0exponentD; isequal(M0,test)


% 2. Convert stressdrop from MPa to bar
sd = round(sd*10);

% 3. Call stf.f
ucmd = sprintf('%sstf sr:%i nn:%i:%i pv:%i rv:%i th:%i k1:%i k2:%i mo:%i:%i sd:%i', ...
    exePath,sr,ns,ne,pv,rv,th,k1,k2,m0mantissaD,m0exponentD,sd);

[~,stdout] = unix(ucmd);

stf      = importdata(outFileFullName);
nsamples = numel(stf); 
t        = (1:nsamples)'*1/sr;