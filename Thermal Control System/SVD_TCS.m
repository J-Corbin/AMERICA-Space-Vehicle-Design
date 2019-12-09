clear, clc

%% Constant Declaration
km_to_m = 1000;
AU_to_m = 1.496e11;
h = 1000 * km_to_m; % Distance to planet, m
Rp= 6378 * km_to_m; % Radius to planet, m
Tp= 250.8; % Blackbody temp of planet, K
rsc= 2; % Radius of spherical spacecraft, m
stefan=5.67e-8; % Stefan-Boltzmann Constant, W/(m^2*K^4)
epsc=0.85; % Average spacecraft emissivity
alsc=0.25; % Average spacecraft absortivity 
SE=1366.1; % Solar constant, W/m^2 
a= 0.3; % Planet albedo
Rsc = (1.496e11 + h) / AU_to_m; % Distance from spacecraft to sun, m
Qint=100;

theta=asin(Rp/(Rp+h)); %%Finding anglular radius of planet as seen by spacecraft, radians
Fscp=0.5*(1-cos(theta)); %% View Factor sc-p
Qsun= alsc*pi*rsc^2*SE*(1/Rsc)^2;
Qalb= alsc*(4*pi*rsc^2)*Fscp*(a*SE*((1/Rsc)^2));
Tsc=(Fscp*Tp^4+((Qsun+Qalb+Qint)/((4*pi*Rsc^2)*stefan*epsc)))^(1/4);

Temp_C = Tsc -273;
fprintf("Spacecraft temperature is:  %f [C]\n", Temp_C)

%% Jack's Junk
%{
theta = asin(Rp/(Rp + h));
rsc = 2;
F = (1 - cos(theta)) / 2;
a = 0.3;
Escp = 0.85;
sigma = 5.67*(10^-8);
alpha = 0.25;
Qi = 10;
AU = 1;
Rsc = 1.496e11;
Se = 1366.1;
Qa = alpha * 4 * pi * rsc^2 * F * (a * Se * (AU/Rsc)^2);
Qs = alpha * pi * rsc^2 * Se * (AU / Rsc)^2;
T = -23 + 273;
Temp = (F * T^4 + (Qs + Qa + Qi)/(4 * pi * rsc^2 * sigma * Escp))^1/4;
%}