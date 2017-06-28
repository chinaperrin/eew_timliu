function [weight] = get_weight(r,m,weightMatrix)

global minimax rEdges mEdges

dr = rEdges(2)-rEdges(1);
dm = mEdges(2)-mEdges(1);

mmin = minimax.mmin;
rmin = minimax.rmin;

ir = floor((r-rmin)/dr)+1;
im = floor((m-mmin)/dm)+1;
    
% Read weight from weight matrix
weight = weightMatrix(ir,im);