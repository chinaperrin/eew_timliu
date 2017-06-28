function [eqs] = get_eqs_subset(eqs,mmin,mmax)

% Return <eqs> structure containing only entries with mmin<=m<mmax

fldNames = fieldnames(eqs);
nfld     = numel(fldNames);
m        = eqs.m;

idx = (m>=mmin & m<mmax);

for ifld = 1:nfld
    eqs.(fldNames{ifld}) = eqs.(fldNames{ifld})(idx);
end