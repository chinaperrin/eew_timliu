function [azimuth, incidence]= azimuth_incidence_ND(Z,E,N,interval)

% Modified again (simplified) by MAM, 170619
%
% Modified script from Nico DeichmaeIdx's incidangle.m
% Finds azimuth and incidence angles from the first motion of an integrated
% seismogram.

% Modified by Falko in 2007 and later by Men-Andrin Meier
% Sept 2013, ETH Zurich

% x,y,z     vertical and horizontal components of displacement time series
% dt        sampling interval (= 1/srate)
% pxIdx     sample-index of p-pick
% tw        time interval over which angles are estimated; time since p-pick

% fit line to horizontal particle motion; assuming x=N,y=E,z=Down (Up?), did not
% test it.
x = N(interval);
y = E(interval);
z = Z(interval);
px = polyfit(x,y,1); % px = polyfit(x(sIdx:eIdx),y(sIdx:eIdx),1);
py = polyfit(y,x,1); % py = polyfit(y(sIdx:eIdx),x(sIdx:eIdx),1);


% calculate azimuth source --> station
angle = atan(mean([1/px(1),py(1)])) * 180/pi;

if (z(end) < 0);              % DOWN
  if (x(end) < 0);
    if (y(end) < 0);
      phi = angle;
    else
      phi = 180 + angle;
    end
  else
    if (y(end) < 0);
      phi = 360 + angle;
    else
      phi = 180 + angle;
    end
  end
else                         % UP
  if (x(end) < 0);
    if (y(end) < 0);
      phi = 180 + angle;
    else
      phi = 360 + angle;
    end
  else
    if (y(end) < 0);
      phi = 180 + angle;
    else
      phi = angle;
    end
  end
end
phi = round(phi);   % azimuth source --> station
azimuth=phi;

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

pr = polyfit(r,z,1); % pr = polyfit(r(1:eIdx),z(1:eIdx),1);
pz = polyfit(z,r,1); % pz = polyfit(z(1:eIdx),r(1:eIdx),1);

% calculate vertical angle of incidence
theta = round(atan(mean([1/pr(1),pz(1)])) * 180/pi);
incidence=theta;















% % plot vertical, radial and transverse displacement seismograms
% if  (o_plot)
%     
%     figure(89); clf;
%     ftSize = 14;
%     ymax = 1.05*max([max(abs(z(sIdx:eIdx))),max(abs(x(sIdx:eIdx))),max(abs(y(sIdx:eIdx)))]);
%     
%     % Map
%     subplot(5,2,[8,10])
%     
%     % Vertical  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
%     subplot(5,2,1:6); hold on
%     plot(t,z,'k')
%     plot(t,x,':r')
%     plot(t,y,'--b')
%     line([t(pxIdx) t(pxIdx)],[-ymax ymax],'color','r','lineWidth',2)
%     l1 = line([t(eIdx) t(eIdx)],  [-ymax ymax],'color',[0 0.4 0],'lineWidth',2);
%     ylabel('[m]','fontSize',ftSize)
%     axis([t(sIdx)-1 t(eIdx)+2 -ymax ymax])
%     legend('Z','N','E')
%     
%    for i = sIdx:eIdx
% 
%        delete(l1)
%        subplot(5,2,[7,9]); hold on; axis equal; hold on
%        l1 = line([t(i) t(i)],[-ymax ymax],'color',[0 0.4 0],'lineWidth',1);
%        plot (x(sIdx:i),y(sIdx:i),'ok')
%        axis([-ymax ymax -ymax ymax])
%        pause(0.01)
%        1+1
%    end
%    
%    plot ([polyval(py,[y(sIdx),y(eIdx)])],[y(sIdx),y(eIdx)],'k')
%    plot ([x(sIdx),x(eIdx)],[polyval(px,[x(sIdx),x(eIdx)])],'k')
% 
%    hold off
%    xlabel ('E-W')
%    ylabel ('N-S')
% end
% 
% if  nargin>4
%    subplot (2,2,3)
%    plot (z,'k')
%    hold on
%    plot (r,'--r')
%    plot (t,':b')
%    hold off
%    %title(evfi)
%    xlabel ('Samples')
%    ylabel ('Displacement')
%    legend ('Z','R','T')
% 
%    % plot vertical vs. radial displacement particle motion
% 
%    subplot (2,2,4)
%    plot (r(1:eIdx),z(1:eIdx),'ok')
%    axis equal
%    hold on
% end
% 
% 
% if  nargin>4
%    plot ([r(1),r(eIdx)],[polyval(pr,[r(1),r(eIdx)])],'k')
%    plot ([polyval(pz,[z(1),z(eIdx)])],[z(1),z(eIdx)],'k')
%    hold off
%    xlabel ('Radial')
%    ylabel ('Vertical')
% end



% write output to file incidangle.log

%fid = fopen('incidangle.log','a');
%fprintf (fid,'%s  %s  %3.0f  %3.0f\n',evfi,datetime,phi,theta);
%fclose(fid);