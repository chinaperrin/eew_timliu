function cmap = shakeMap_cscale(nc)

clrMAP = 1/256*[255	255	255;
    255	255	255; ...
    191	204	255; ...
    160	230	255; ...
    128	255	255; ...
    122	255	147; ...
    255	255	0; ...
    255	200	0; ...
    255	145	0; ...
    255	0	0; ...
    200	0	0];

intVect_orig   = (0:10)';
intVect_target = linspace(intVect_orig(1),intVect_orig(end),nc);

cmap      = zeros(nc,3);
cmap(:,1) = interp1(intVect_orig,clrMAP(:,1),intVect_target);
cmap(:,2) = interp1(intVect_orig,clrMAP(:,2),intVect_target);
cmap(:,3) = interp1(intVect_orig,clrMAP(:,3),intVect_target);