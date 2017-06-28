function [arrayOut] = double2single_cellArray(arrayIn)
%arrayIn=ssrte.pythia2_ss.immi.prct;

arrayOut = cell(size(arrayIn));
for ii = 1:numel(arrayIn)
    arrayOut{ii} = single(arrayIn{ii});
end
%whos arrayIn arrayOut
%a = arrayOut{12};
%b = arrayIn{12};
%whos a b
%save('arrayIn','arrayIn')
%save('arrayOut','arrayOut')