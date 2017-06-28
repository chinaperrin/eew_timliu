function [azimuth, incidence]= azimuth_incidence(x,y,z,dt,plot_option)

% this is a modified versoin of Nico Deichmann's incidangle.m
% the unix part was removed, input is now 3 components of a seismogram
% if no zero crossing is found the steepest slope of the seismogram is
% taken instead (max value acceleorgram)
% a filtered velocity seismogram as input is recommended
%
% Falko Bethmann
% ETH Zurich
% Jan 2006
% plot option changed to nargin
% May 2007
% mean substracted for integration
% Aug 2007

if size(x,2)==1
x = x';
y = y';
z = z';
end


%if max(z)<0 || min(z)>0
    %no zero crossing
    z=z-mean(z);
%end
%if max(y)<0 || min(y)>0
    %no zero crossing
    y=y-mean(y);
%end
%if max(x)<0 || min(x)>0
    %no zero crossing
    x=x-mean(x);
%end
% % integrate to displacement only in case of accelerometer data
% 
% x=bandpass(x,1,10,dt,2)
% y=bandpass(y,1,10,dt,2)
% z=bandpass(z,1,10,dt,2)
% 
% z=z-mean(z);
% y=y-mean(y);
% x=x-mean(x);
% 
% z = cumtrapz(z) .* dt;
% y = cumtrapz(y) .* dt;
% x = cumtrapz(x) .* dt;


% plot velocity seismograms
if  nargin>4
   figure(989); clf;
   subplot (2,2,1)
   plot(z,'k')
   hold on
   plot(y,'--r')
   plot(x,':b')
   hold off
   %title (datetime)
   xlabel ('Samples')
   %ylabel (units)
   legend ('Z','N','E')   
end

% find first zero-crossing of velocity seismograms
% if no crossing found take max value of accelerogram
acc_z=gradient(z,dt);
acc_x=gradient(x,dt);
acc_y=gradient(y,dt);




n = length(z);
i  = find(z(1:n-1).*z(2:n) < 0);
    if isempty(i)
        [min_z,iz]=min(acc_z);
    else
        iz = i(1);   
    end

if (iz < 2)
    iz = i(2);
end

i  = find(y(1:n-1).*y(2:n) < 0);

    if isempty(i)
        [min_y,iy]=min(acc_y);
    else
        iy = i(1);   
    end
    
if (iy < 2)
    iy = i(2);
end

i  = find(x(1:n-1).*x(2:n) < 0);
    if isempty(i)
        [min_x,ix]=min(acc_x);
    else
        ix = i(1);   
    end
if (ix < 2)
    ix = i(2);
end


nn = min([iz,iy,ix]) + 1;

% integrate to displacement

z=z-mean(z);
y=y-mean(y);
x=x-mean(x);

z = cumtrapz(z(1:n)-z(1,1)) .* dt;
y = cumtrapz(y(1:n)-y(1,1)) .* dt;
x = cumtrapz(x(1:n)-x(1,1)) .* dt;





% plot NS vs. EW displacement particle motion
if  nargin>4
   subplot (2,2,2)
   plot (x(1:nn),y(1:nn),'ok')
   axis equal
   hold on
end

% fit line to horizontal particle motion

px = polyfit(x(1:nn),y(1:nn),1);
py = polyfit(y(1:nn),x(1:nn),1);

if  nargin>4
   plot ([polyval(py,[y(1),y(nn)])],[y(1),y(nn)],'k')
   plot ([x(1),x(nn)],[polyval(px,[x(1),x(nn)])],'k')

   hold off
   xlabel ('E-W')
   ylabel ('N-S')
end

% calculate azimuth source --> station

angle = atan(mean([1/px(1),py(1)])) * 180/pi;

if (z(nn) < 0);              % DOWN
  if (x(nn) < 0);
    if (y(nn) < 0);
      phi = angle;
    else
      phi = 180 + angle;
    end
  else
    if (y(nn) < 0);
      phi = 360 + angle;
    else
      phi = 180 + angle;
    end
  end
else                         % UP
  if (x(nn) < 0);
    if (y(nn) < 0);
      phi = 180 + angle;
    else
      phi = 360 + angle;
    end
  else
    if (y(nn) < 0);
      phi = 180 + angle;
    else
      phi = angle;
    end
  end
end
phi = round(phi);   % azimuth source --> station
azimuth=phi;

if  nargin>4
   title (['phi = ',num2str(phi)])
end

% rotate horizontal components to radial and transverse

sinphi = sin(phi * .0174533);
cosphi = cos(phi * .0174533);
a(1,1) =  sinphi;
a(1,2) =  cosphi;
a(2,1) =  cosphi;
a(2,2) = -sinphi;

r = a(1,:)*[x;y];
r = r;
t = a(2,:)*[x;y];
t = t;

% plot vertical, radial and transverse displacement seismograms

if  nargin>4
   subplot (2,2,3)
   plot (z,'k')
   hold on
   plot (r,'--r')
   plot (t,':b')
   hold off
   %title(evfi)
   xlabel ('Samples')
   ylabel ('Displacement')
   legend ('Z','R','T')

   % plot vertical vs. radial displacement particle motion

   subplot (2,2,4)
   plot (r(1:nn),z(1:nn),'ok')
   axis equal
   hold on
end

% fit line to vertical particle motion

pr = polyfit(r(1:nn),z(1:nn),1);
pz = polyfit(z(1:nn),r(1:nn),1);

if  nargin>4
   plot ([r(1),r(nn)],[polyval(pr,[r(1),r(nn)])],'k')
   plot ([polyval(pz,[z(1),z(nn)])],[z(1),z(nn)],'k')
   hold off
   xlabel ('Radial')
   ylabel ('Vertical')
end

% calculate vertical angle of incidence

theta = round(atan(mean([1/pr(1),pz(1)])) * 180/pi);

if  nargin>4
   title (['theta = ',num2str(theta)])
end
incidence=theta;

% write output to file incidangle.log

%fid = fopen('incidangle.log','a');
%fprintf (fid,'%s  %s  %3.0f  %3.0f\n',evfi,datetime,phi,theta);
%fclose(fid);
