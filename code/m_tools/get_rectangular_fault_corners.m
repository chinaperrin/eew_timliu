function plane = get_rectangular_fault_corners(clat,clon,cz,strike,dip,L,W,fig)
% Computes 3D coordinates of the four corners of a square fault for a 
% specified reference corner point at (clat,clon,cz) 
%
% INPUT UNITS: all angles [deg], L,W & z [km]
%
% menandrin@gmail.com, 160929

% Conversions ------------------------
km2lat  = 1/111.2;           % km2deg
r_earth = 6371.009;
lon2km  = r_earth*cos(clat*pi/180)*pi/180;
km2lon  = 1/lon2km;
strike  = strike*pi/180;     % deg2rad
dip     = dip   *pi/180;
% ____________________________________

nl = numel(L);
nw = numel(W);
LL = repmat(L,nw,1);
WW = repmat(W,1,nl);

p.x = zeros(nw,nl);
p.y = zeros(nw,nl);
p.z = zeros(nw,nl);

dx1 = LL*sin(strike);
dy1 = LL*cos(strike);
dz  = WW*sin(dip);
dh  = WW*cos(dip);
dx2 = dh*cos(strike);
dy2 = dh*sin(strike+pi);

plane.x   = 0+dx1+dx2;                 
plane.y   = 0+dy1+dy2;
plane.z   = cz+dz;
plane.lat = clat+(dy1+dy2)*km2lat;
plane.lon = clon+(dx1+dx2)*km2lon;

p1.x = plane.x(1  ,  1); p2.x = plane.x(1  ,end); p3.x = plane.x(end,  1); p4.x = plane.x(end,end);
p1.y = plane.y(1  ,  1); p2.y = plane.y(1  ,end); p3.y = plane.y(end,  1); p4.y = plane.y(end,end);
p1.z = plane.z(1  ,  1); p2.z = plane.z(1  ,end); p3.z = plane.z(end,  1); p4.z = plane.z(end,end);
p1.lat = plane.lat(1  ,  1); p2.lat = plane.lat(1  ,end); p3.lat = plane.lat(end,  1); p4.lat = plane.lat(end,end);
p1.lon = plane.lon(1  ,  1); p2.lon = plane.lon(1  ,end); p3.lon = plane.lon(end,  1); p4.lon = plane.lon(end,end);

if fig.plot
    xp = [p1.x p2.x p4.x p3.x];
    yp = [p1.y p2.y p4.y p3.y];
    zp = [p1.z p2.z p4.z p3.z];
    
    lonp = [p1.lon p2.lon p4.lon p3.lon];
    latp = [p1.lat p2.lat p4.lat p3.lat];
    
    figure(12); clf;
%     subplot(2,1,1); hold on; grid on; box on;
%     patch(xp,yp,'b')
%     fill3(xp,yp,zp,'r')
%     set(gca,'xlim',[-1.2*L(end) 1.2*L(end)],'ylim',[-1.2*L(end) 1.2*L(end)],'zlim',[0 1.2*(cz+max(dz(:)))],'zdir','reverse')
%     view(fig.view)
%     xlabel('x'); ylabel('y'); zlabel('z')
%     plot3([0 0], get(gca,'ylim'),[0 0],'-k')
%     plot3(get(gca,'xlim'), [0 0],[0 0],'-k')
%     
%     subplot(2,1,2); 
    hold on; grid on; box on;
    patch(lonp,latp,'b','faceAlpha',0.3)
    fill3(lonp,latp,zp,'r')
    Ldeg = L*km2lat; 
    set(gca,'xlim',[-2*L(end) 2*L(end)],'ylim',[-2*L(end) 2*L(end)],'zlim',[0 2*(cz+max(dz(:)))],'zdir','reverse')
    set(gca,'xlim',[clon-2*Ldeg(end) clon+2*Ldeg(end)],'ylim',[clat-2*Ldeg(end) clat+2*Ldeg(end)],'zlim',[0 2*(cz+max(dz(:)))],'zdir','reverse')
    view(fig.view)
    xlabel('lon'); ylabel('lat'); zlabel('z')
    plot3([clon clon], get(gca,'ylim'),[0 0],'-k')
    plot3(get(gca,'xlim'), [clat clat],[0 0],'-k')    
end