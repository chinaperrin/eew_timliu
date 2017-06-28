function [nWorkers] = startPool(nWorkersSuggestion)

uname = getenv('USER');
if strcmp(uname,'mameier'); nWorkers = 2; else nWorkers = nWorkersSuggestion; end

if matlabpool('size')==0            % checking to see if my pool is already open
    matlabpool('open',nWorkers);
else
    fprintf(1,['Matlabpool already open (',num2str(matlabpool('size')),' workers)\n'])
end
