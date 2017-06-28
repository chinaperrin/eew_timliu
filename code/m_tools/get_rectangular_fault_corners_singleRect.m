function rect = get_rectangular_fault_corners_singleRect(clat,clon,cz,strike,dip,L,W,fig)
% Computes 3D coordinates of the four corners of a square fault for a 
% specified reference corner point at (clat,clon,cz) 
%
% INPUT UNITS: all angles [deg], L,W & z [km]
%
% menandrin@gmail.com, 160929


% Conversions -----------------------
km2lat = 1/111.2;           % km2deg
km2lon = 1/lon2km(1,clat);
strike = strike*pi/180;     % deg2rad
dip    = dip   *pi/180;
% ___________________________________

dx1    = L*sin(strike);
dy1    = L*cos(strike);
dh     = W*cos(dip);
dx2    = dh*cos(strike);
dy2    = dh*sin(strike+pi);
dz     = W*sin(dip);

p1.x   = 0;                 % Reference point
p1.y   = 0;
p1.z   = cz;
p1.lat = clat;
p1.lon = clon;

p2.x   = p1.x+dx1;          % Second point at surface
p2.y   = p1.y+dy1;
p2.z   = cz;
p2.lon = clon+dx1*km2lon;
p2.lat = clat+dy1*km2lat;

p3.x   = p1.x+dx2;          % Point at depth beneath reference point
p3.y   = p1.y+dy2;
p3.z   = cz+dz;
p3.lon = clon+dx2*km2lon;
p3.lat = clat+dy2*km2lat;

p4.x   = p1.x+dx1+dx2;      % Second point at depth
p4.y   = p1.y+dy1+dy2;
p4.z   = cz+dz;
p4.lon = clon+(dx1+dx2)*km2lon;
p4.lat = clat+(dy1+dy2)*km2lat;

rect.p1 = p1;
rect.p2 = p2;
rect.p3 = p3;
rect.p4 = p4;



if fig.plot
    xp = [p1.x p2.x p4.x p3.x];
    yp = [p1.y p2.y p4.y p3.y];
    zp = [p1.z p2.z p4.z p3.z];
    
    lonp = [p1.lon p2.lon p4.lon p3.lon];
    latp = [p1.lat p2.lat p4.lat p3.lat];
    
    figure(12); clf;
    subplot(2,1,1); hold on; grid on; box on;
    patch(xp,yp,'b')
    fill3(xp,yp,zp,'r')
    set(gca,'xlim',[-1.2*L 1.2*L],'ylim',[-1.2*L 1.2*L],'zlim',[0 1.2*(cz+dz)],'zdir','reverse')
    view(fig.view)
    xlabel('x'); ylabel('y'); zlabel('z')
    plot3([0 0], get(gca,'ylim'),[0 0],'-k')
    plot3(get(gca,'xlim'), [0 0],[0 0],'-k')
    
    subplot(2,1,2); hold on; grid on; box on;
    patch(lonp,latp,'b')
    fill3(lonp,latp,zp,'r')
    Ldeg = L*km2lat; 
    set(gca,'xlim',[clon-1.2*Ldeg clon+1.2*Ldeg],'ylim',[clat-1.2*Ldeg clat+1.2*Ldeg],'zlim',[0 1.2*(cz+dz)],'zdir','reverse')
    view(fig.view)
    xlabel('lon'); ylabel('lat'); zlabel('z')
    plot3([clon clon], get(gca,'ylim'),[0 0],'-k')
    plot3(get(gca,'xlim'), [clat clat],[0 0],'-k')
    
end