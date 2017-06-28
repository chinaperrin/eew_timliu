function Td = get_srcDuration_HanksThatcher72(Mw)

% ?? = 8.47*M0 (fc/beta)^3  % Hanks and Thatcher (1972), Baltay and Hanks (2015) eq. 5
% stress drop = 4.64 MPa bottom of p. 2856 (Model Development)

beta        = 3.5*1e3; % km/s -> m/s
delta_sigma = 5*1e6;   % MPa -> Pa = N/m^2

%Mw=2:.5:8; %Mw=linspace(4.5,9,30);
%Mo=Mw2Mo(Mw)*1e-7; % dyne*cm -> N*m
Mo = magnitude2moment(Mw);

fc = beta * [delta_sigma./(8.47*Mo)].^(1/3);
Td = 1./fc;