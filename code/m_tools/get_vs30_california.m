function [vs30,drVect] = get_vs30_california(stationList)
% READ VS30 VALUE FOR LIST OF STATIONS FROM NEAREST GRID POINT ON VS30-MAP
%

socal_vs30_fullName = '~/Documents/eq/data/socal/vs30/usgs_cgs_geology_60s.txt';
tmp                 = importdata(socal_vs30_fullName);

stLonVect = tmp(:,1); 
stLatVect = tmp(:,2); 
vs30Vect  = tmp(:,3); 
n         = numel(stLatVect);

% Plot and print map
% addpath(genpath('../../../matlab/colorscales/'))
% clf; scatter(stLonVect,stLatVect,15,vs30Vect,'filled')
% clrIncs  = 170:10:800;
% clrSteps = [0;180;240;300;360;490;620;760;800];
% clrmap   = make_cScale_w_irregular_steps(clrIncs,clrSteps);
% colorbar; colormap(clrmap)
% set(gca,'xlim',[-125 -114],'ylim',[31.5 42])
% print('-dpng','~/Documents/eq/data/socal/vs30/usgs_cgs_geology_60s.txt.png')


nst    = numel(stationList.lat);
vs30   = zeros(nst,1);
drVect = zeros(nst,1);
for ist = 1:nst
    
    print_iteration_numbers(ist,nst,'hundreds')
    
    stLat = stationList.lat(ist);
    stLon = stationList.lon(ist);

    dlat = stLat-stLatVect;
    dx   = deg2km(dlat);
    dlon = stLon-stLonVect;
    dy   = lon2km(dlon,stLat*ones(n,1));
    dr   = sqrt(dx.^2 + dy.^2);
    
    [valMin,idxMin] = min(dr);
    vs30(ist)       = vs30Vect(idxMin);
    drVect(ist)     = valMin;
end
fprintf(1,sprintf('VS30-map: Max distance to closest grid point was %5.1fkm',max(drVect)))















% _________________________________________________________________________
function [lon_km] = lon2km(lon_deg,lat_deg)
% Converts relative longitude values <lon_deg> into kilometer as a function 
% of the respective latitudes <lat_deg>. Assumes a sphere.
% mam, 120427

r_earth = 6371.009;

if ( ~isequal(size(lon_deg),size(lat_deg)))
    error('<lon_deg> and <lat_deg> must have the same format')
end

deg2km_lon = r_earth*cos(lat_deg*pi/180)*pi/180;
lon_km     = lon_deg.*deg2km_lon;
end




% _________________________________________________________________________
function print_iteration_numbers(itr,ntr,printLevel)

%eins      = linspace(1e0,1e4,1e4);
tens      = linspace(1e1,1e4,1e3);
hundreds  = linspace(1e2,1e5,1e3);
thousands = linspace(1e3,1e6,1e3);

% EXAMPLE
% ni = 120;
% for ii=1:ni
%     print_iteration_numbers(ii,ni,'ones')
% end

%if nargin<4; o_printTime = false; end
if itr==1; tic; end


% Print hundreds on same line and start new line every 1,000th iteration 
if strcmp(printLevel,'hundreds')
    if ismember(itr,[1,hundreds,ntr]);
        %if ismember(itr,[1,thousands]); fprintf(1,sprintf('%16.4fs\n%i -- ',toc,ntr)); end
        if ismember(itr,[1,thousands]); fprintf(1,sprintf('\n%i -- ',ntr)); end
        fprintf(1,sprintf('%i .. ',itr));
        if itr==ntr; fprintf(1,' done.\n'); end
    end
    
    
% Print tens on same line and start new line every 100th iteration
elseif strcmp(printLevel,'tens')
    if ismember(itr,[1,tens,ntr]);
        %if ismember(itr,[1,hundreds]);  fprintf(1,sprintf('%16.4fs\n%i -- ',toc,ntr)); end
        if ismember(itr,[1,hundreds]); fprintf(1,sprintf('\n%i -- ',ntr)); end
        fprintf(1,sprintf('%i .. ',itr));
        if itr==ntr; fprintf(1,' done.\n'); end
    end
    
% Print singles on same line and start new line every 10th iteration 
elseif strcmp(printLevel,'ones')
        %if ismember(itr,[1,tens]); fprintf(1,sprintf('%16.4fs\n%i -- ',toc,ntr)); end
        if ismember(itr,[1,tens]); fprintf(1,sprintf('\n%i -- ',ntr)); end
        fprintf(1,sprintf('%i .. ',itr));
    end
end
end
