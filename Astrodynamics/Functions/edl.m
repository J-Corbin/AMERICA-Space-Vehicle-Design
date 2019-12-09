clear, clc
global beta r u hscale rho0 mass

%% Calculating Ballistic Coefficient
% Constants ballistic coefficient (beta) depends on
m = 1; %kg
cd = 1.6; %drag coefficient, taken from MER capsule
A = 1; %reference area (m^2)

beta = m / (cd * A); %Ballistic coefficient
%% Allen-Eggers for Ballistic Entry
%{
% Constants declaration
gamma = 1.3; %ratio of specific heats for Mars
H = 1; %atmospheric scale height for Mars

C = - p_0 * H / (2 * beta * sin(gamma));
v = v_atm * exp(C * exp(-h / H));
%}
r = 3376.2; %Mars' radius
u = 42828; %Mars' gravitational parameter
hscale = 11.1; %Mars' scale height
mass = 700; %mass of entry vehicle
rho0 = .02; %Mars' surface density

tspan = [0 inf];
f0 = [4.0, r * 2, 10 * pi / 180];
Opt = odeset('Events', @landed);
[time, f0] = ode45(@ballistic, tspan, f0, Opt);

%% Numerical Solution to Ballistic Entry Equations of Motion
function dydt = ballistic(t, y)
    global r u hscale rho0 mass beta

    dydt = zeros(size(y));
    
    v = y(1);
    h = y(2);
    gamma = y(3);
    
    % Calculate gravitational acceleration
    g = -(mass * u) / (r + h)^2;
    
    % Calculate density
    rho = rho0 * exp(-h / hscale);
    
    % Differential equations
    v_dot = -rho * v^2 / (2 * beta) + g * cos(gamma);
    h_dot = -v * sin(gamma);
    gamma_dot = (-v^2 * cos(gamma) / (r + h) + g * cos(gamma)) / v;
    
    %{
    %Intermediate calculations
    gravity = -(mass * u) / (r_earth + h)^2;
    density = rho0 * exp(-h / hscale);
    drag = -0.5 * density * Cd * A * v^2;
    
    %Differential equations
    phi_dot = (u / (r_earth + h)^2) * sin(psi) / v;
    v_dot = (T + drag + gravity * cos(psi)) / mass;
    h_dot = v * cos(psi);
    x_dot = v * sin(psi);
    theta_dot = (v * sin(psi)) / (r_earth + h);
    psi_dot = phi_dot - theta_dot;
    %}
    
    dydt(1) = v_dot;
    dydt(2) = h_dot;
    dydt(3) = gamma_dot;
end

function [value,isterminal,direction] = landed(~,y)
global r
%Event funtion to stop integration when lander reaches the surface
if y(2) > r
 value = 1; %Keep going
else %If not
 value = 0; %Then stop
end
isterminal = 1; %Terminate integration when condtion met
direction = 0; %Direction doesn't matter
end