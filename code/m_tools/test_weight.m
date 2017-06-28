clear all

global mmin mmax dmin dmax

mmin     = 1;
mmax     = 7;
dmin     = 0;
dmax     = 50;


% data
M  = [1.4, 3.2, 3.8, 5.4, 5.3, 5.5,5.8];
D  = 10*[3.1, 2.1,2.9, 1.3, 1.3, .17, 1.2];


% compute weight matrix
[weightmat,wm,wd] = get_weight_matrix(M,D,1);

% read weight of particular datum
m = 1.5;
d = 30;

[weight] = get_weight(m,d,weightmat,wm,wd);


%histmat = hist2(xdata,ydata,xedges,yedges)