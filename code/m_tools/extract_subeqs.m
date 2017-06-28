function eqsnew = extract_subeqs(eqs,mminTarget);

keepMe = eqs.m>mminTarget;

fdnameList = fieldnames(eqs);
nfields    = numel(fdnameList);
for ifield = 1:nfields
    fdName = fdnameList{ifield};
    eqsnew.(fdName) = eqs.(fdName)(keepMe);
end