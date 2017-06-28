function out = sampleVect2pdf(aVal,aVarIn,pctLevels)

%% Make histogram of GMPs, normalise to pdf
if numel(aVarIn)==1;  da     = aVarIn;
                      p01    = prctile(aVal,.1);
                      p99    = prctile(aVal,99.9);
                      aEdges = p01:da:p99;
else numel(aVarIn)>1; aEdges = aVarIn;
                      p01    = prctile(aVal,.1);    % dummies for output
                      p99    = prctile(aVal,99.9);  % dummies for output
end
na   = numel(aEdges);
nval = numel(aVal);

if na>1 %if p01~=p99
    
    nc = zeros(na-1,1);
    for ia = 1:na-1
        aLo    = aEdges(ia);
        aUp    = aEdges(ia+1);
        nc(ia) = sum(aVal>=aLo &aVal<aUp);
    end
    %pA     = nc./sum(nc);
    pA     = nc./nval;
    
    %% Compute maximum probability estimate
    [pmax,imax] = max(pA);
    mpe         = 1/2*(aEdges(imax)+aEdges(imax+1));
    
else
    pA     = 1;
    aEdges = [p01-da/2,p01+da/2];
    mpe    = p01;
    pmax   = 1;
    nc     = numel(aVal);
end

% Optionally ompute percentiles
if nargin>2
    out.pct = prctile(aVal,pctLevels);
end

% Prepare output
out.p      = pA;
out.counts = nc;
out.aEdges = aEdges;
out.mpe    = mpe;
out.pmax   = pmax;
out.nval   = nval;
out.fOutLo = sum(aVal<p01)/nval;
out.fOutUp = sum(aVal>=p99)/nval;