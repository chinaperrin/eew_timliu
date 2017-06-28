%% OVERVIEW
%
% How to estimate durations
% 1. Olson and Allen, 2005, Nature; rupture length from scaling laws,
%    assume rupture velocity --> estimate duration
% 2. Regression of source durations from STF catalogs
% 3. Scaling law between moment and duration, calibrated at some magnitude $
%    for which we know duration
% 4. Scaling law between length/vr and duration, use Wells and Coppersmith,
%    1994, BSSA
% 5. From theoretical source model, e.g. Sato / Hirasawa

M  = (2:.1:9)';
nm = numel(M);
ftSize = 15;

figure(821); clf; hold on; grid on
xlabel('M_w','fontSize',ftSize)
ylabel('Source Duration [sec]','fontSize',ftSize)


% 1. ----------------------------------------------------------------------
% Olson and Allen, 2005, Nature Only vaguely describe how they compute 
% source duration. From their Figure 
% M3.0: duration D = 0.1sec
% M5.2: duration D = 2sec
% log10(D) = a+b*M

% Manual
b    = ( log10(.1) - log10(2)) / -2.2;
a    = 1/2 * ( log10(.1) + log10(2) -8.2*b);
D3p0 = 10^(a + b*3.0);
D5p2 = 10^(a + b*5.2);
D1   = 10.^(a + b*M);

% Inverse
C    = [1 1; 3 5.2];
invC = inv(C);
ab   = [log10(.1); log10(2)]'*invC;
a2   = ab(1); 
b2   = ab(2);

gcf; 
p1 = plot(M,D1,'lineWidth',2,'color','k');
set(gca,'yscale','log','fontSize',ftSize)



% 2. ----------------------------------------------------------------------

% 
% % 3. ----------------------------------------------------------------------
% vr         = 3000;           % rupture velocity in [m/s]
% stressdrop = 1e6;            % Stress drop in [Pa];
% cc         = 1.5*1e3;           % Constant
% c          = 1/(vr*stressdrop^(1/3));
% D2         = cc*c*10.^(.5*M);
% gcf; 
% p2 = plot(M,D2,'lineWidth',2,'color','r');
% 



% 4. ----------------------------------------------------------------------
% Wells and Coppersmith regressions for subsurface rupture length and
% magnitude; using all mechanisms.
%log10(RLD) = a + b*M
vr = 2700;           % rupture velocity in [m/s]
a  = -2.44;
sa = 0.11;
b  = 0.59;
sb = 0.02;
RLD = 10.^(a + b*M)*1e3;        % Subsurface rupture length in [m]
D3  = RLD./vr;
gcf; 
p3 = plot(M,D3,'lineWidth',2,'color','m');



% 4b. Get duration distribution by sampling a, b, vr 
n   = 1e5;
ar  = a + sa*randn(n,1);
br  = b + sb*randn(n,1);
vrr = 2.4 + 0.6*rand(n,1);
MM  = repmat(M',n,1);
bbr = repmat(br,1,nm);
aar = repmat(ar,1,nm);
vvrr = repmat(vrr,1,nm); 
logRLD = (aar + bbr.*MM);
RLD    = 10.^logRLD;
D3r = RLD./vvrr;
pct = prctile(D3r,[16 50 84]);

gcf;
p4 = plot(M,pct(2,:),'r','lineWidth',2);
plot(M,pct(1,:),'r','lineWidth',1)
plot(M,pct(3,:),'r','lineWidth',1)



% 5. ----------------------------------------------------------------------


% Finish plot
l1 = legend([p1; p4],'Olson and Allen 2005','D ~ L/vr; Wells and Coppersmith, 1994');
set(l1,'location','northwest','fontSize',ftSize)