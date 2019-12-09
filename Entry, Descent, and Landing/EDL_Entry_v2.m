%% EDL Entry Code
% Jasmina Muftic
% Space Vehicle Design Capstone
% Fall 2019

%% Constants
v_atm    = 5.94*1000;       % Entry velocity [m/s]
gam0     = deg2rad(13);     % Flight Path Angle [rad]
m        = 700;             % Entry mass [kg]
Cd       = 1.65;            % Drag coefficient of 70 deg sphere-cone heatshield
g0_mars  = 3.711;           % Mars surface gravity [m/s^2]
g0_earth = 9.81;            % Earth surface gravity [m/s^2]
T_0      = 240;             % Estimate for average martian surface temperature [K]
R        = 192;             % Specific Martian gas constant [J/kg*K]
radius   = 2;               % Cross-sectional radius
A        = pi*radius^2;     % Cross-sectional area [m^2]
r_p      = 3389.5*1000;     % Martian radius [m]
rho_0    = 0.02;            % Martian surface density [kg/m^3]
spec_h   = 1.3;             % Ratio of specific heats
Cx       = 1.3;             % Opening force coefficient

%% Entry Vehicle Model Variables
h_0      = 125*1000;        % Initial altitude [m]
beta     = m./(Cd.*A);      % Ballistic coefficient
H_mars   = 11100;           %(R*T_0/g0_mars); % Atmospheric scale height [m]

f = @(t,x) [(-(rho_0 * exp(-x(2)/H_mars)) * x(1)^2)./(2.*beta) + g0_mars * sin(x(3)); -x(1)*sin(x(3)); (-x(1)*cos(x(3)))/(r_p + x(2)) + (g0_mars * cos(x(3)))/x(1)];

[t1, x1] = ode45(f, [0 370], [v_atm h_0 gam0]);

v1 = x1(:,1);                   % descent velocity [m/s]
descent  = x1(:,2);             % descent height [m]
fpa      = rad2deg(x1(:,3));    % flight path angle [degrees]

C          = (-rho_0 * H_mars)/(2*beta*sin(gam0));
peak_decel = (v_atm^2 * sin(gam0))/(53.3*H_mars); % Peak deceleration in mars g's[m/s^2]
alt_decel  = H_mars*log(-2*C);                    % Altitude of peak deceleration [m]
dydx       = gradient(v1(:)) ./ gradient(t1(:));  % Acceleration
horiz_d    = cumtrapz(t1,(v1.*cosd(fpa)));        % Horizontal distance traveled

%% Parachute Deployment
Cd2       = 0.5;             % Drag coefficient of 70 deg sphere-cone heatshield
radius2   = 12.8/2;          % Cross-sectional radius
A2        = pi*radius2^2;    % Cross-sectional area [m^2]
p_deploy = 33;
t_interv = 300;

tspan = [t1(end-p_deploy), t1(end-p_deploy)+t_interv];
fpa0 = deg2rad(fpa(end-p_deploy));
x0 = v1(end-p_deploy)*sin(fpa0);
x0_1 = v1(end-p_deploy)*cos(fpa0);

[t2,x2]      = ode45(@(t2,x2) (m*g0_mars - 0.5*Cd2*rho_0*A2*x2^2)/m, tspan, x0); % Vertical component of velocity
[t3,x3]      = ode45(@(t3,x3) (-0.5*Cd2*rho_0*A2*x3^2)/m, tspan, x0_1);          % Horizontal component of velocity

addon1 = zeros(abs((length(x2)-length(x3))),1);
addon2 = x2(end) .* ones(abs((length(x2)-length(x3))),1);

if length(x2)>length(x3)
    x3 = [x3; addon1];
else
    x2 = [x2; addon2];
    t2 = t3;
end

fpa2 = atan(x2./x3);                                % Flight path angle
v2 = sqrt(x2.^2 + x3.^2);                           % Velocity (magnitude)
dydx2 = gradient(v2(:)) ./ gradient(t2(:));         % Acceleration
descent2 = descent(end-p_deploy) - cumtrapz(t2,x2); % Altitude
horiz_d2 = cumtrapz(t2,x3) + horiz_d(end-p_deploy); % Horizontal distance traveled

terminal = sqrt((2*m*g0_mars)/(Cd2*rho_0*A2));      % Expected terminal descent

% Combining pre- and post- deployment parameters
t = vertcat(t1(1:end-p_deploy),t2(2:end));              % Time [s]
v = vertcat(v1(1:end-p_deploy),v2(2:end));              % Velocity [m/s]
a = vertcat(dydx(1:end-p_deploy),dydx2(2:end));         % Acceleration [m/s^2]
d = vertcat(descent(1:end-p_deploy),descent2(2:end));   % Altitude [m]
y = vertcat(fpa(1:end-p_deploy),rad2deg(fpa2(2:end)));  % Flight path angle [deg]
x = vertcat(horiz_d(1:end-p_deploy),horiz_d2(2:end));   % Horizontal distance [m]

pressure = zeros(numel(d),1);
T = zeros(numel(d),1);
for i = 1:1:numel(d)
        T(i) = (-31 - (0.000998.*d(i))) + 273; % Atmospheric Temperature [K]
end 

for i = 1:1:numel(d)
        pressure(i) = 0.699 * exp(-0.00009*d(i)); % Pressure [Pa]
end

rho = pressure./(0.1921 .* T); % Density 
max_load = 0.5*rho(length(v1)-p_deploy)*(v1(end-p_deploy))^2*Cd2*A2*Cx; % Peak Loading [N]

%% Plots

plot(t,v)
xlabel('Time [s]')
ylabel('Velocity [m/s]')

figure()
plot(t, a)
xlabel('Time [s]')
ylabel('Deceleration [m/s^2]')

figure()
plot(t, d./1000)
xlabel('Time [s]')
ylabel('Height [km]')
grid on

figure()
plot(t, y)
xlabel('Time [s]')
ylabel('Flight Path angle [degrees]')

figure()
plot(v./1000,d./1000)
xlabel('Velocity [km/s]')
ylabel('Altitude [km]')
grid on

figure()
plot(t,x./1000);
xlabel('Time [s]')
ylabel('Lateral Distance Traveled [km]')


