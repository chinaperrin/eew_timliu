function [bAz,ray] = get_backAzimuth(eqLat,eqLon,stLat,stLon)
% Compute the back-azimuth (bAz) and store a ray-coordinates for plotting
% Unit: [deg]

dlat = stLat-eqLat;
dlon = stLon-eqLon;

fprintf(1,'not sure if this script returns the correct azimuth. Check before using.\n')
% assuming a plane ...
alpha = rad2deg(atan(dlon/dlat));
if (dlon >= 0)
    az = 90 - alpha;
elseif (dlon < 0)
    az = 270 - alpha;
end

% az is the azimuth from epiC to station, now get BACK-azimuth from station to epiC:
bAz = az + 180;
if (bAz >= 360); bAz = bAz-360; end

% Compute rays for plotting
if (dlat > 0)
    ray.lat = linspace(stLat,stLat - 10,100);
else
    ray.lat = linspace(stLat,stLat + 10,100);
end
ray.lon = (ray.lat - stLat).*tan(deg2rad(bAz)) + stLon;
