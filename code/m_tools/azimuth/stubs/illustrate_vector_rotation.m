v = [-1 -12];
dtheta = 5;
theta  = 0;
n      = 360/dtheta;

figure(22); clf; hold on
quiver(0,0,v(1),v(2))


for i = 1:n

    theta = theta+dtheta;

    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    vR = v*R;

    quiver(0,0,vR(1),vR(2))
end