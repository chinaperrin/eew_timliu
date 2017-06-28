function Dp=bt2015_path_duration(R)
% Dp(sec)=bt2015_path_duration(Rps_km)
%
% Path duration from Boore and Thompson (2015)
% Duration is in seconds.
% Rps is point-source distance in km
% These durations are for active crustal regions (ACRs) (Table 2).

TABLE=...
[0.0   0.0
7.0   2.4
45.0  8.4
125.0 10.9
175.0 17.4
270.0 34.2];

x=TABLE(:,1);
y=TABLE(:,2);

Dp=interp1(x,y,R,'linear',nan);

ii=find(R>x(end));
Dp(ii)=y(end) + 0.156*(R(ii)-x(end));

return
