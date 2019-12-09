clear, clc 

%% Setting constants
global g

% Set gravity acceleration value
g = -3.711; % gravity acceleration for Mars

%% Setting initial conditions and timespan for integration
% Set arbitrarily long timespan to integrate over
tspan=[0, 10000];

% Create initial condition array
r = [100, 100, 100]; % Initial position coordinates (m)
v = [-0.1, -0.1, -0.1]; % Initial velocity values (m/s)
m = 435; % Initial mass (kg)

IC = [r, v, m]; % Array of initial conditions to be fed into ode45

%% Perform ode45 integration
% Specify optional event to stop integration once spacecraft has reached
% the ground (z = 0)
Opt1 = odeset('Events', @Ground);

% Call ode45 with the initial conditions, timespan, and event for given
% function
[t_out,state_out]=ode45(@odefun,tspan,IC, Opt1);


%% Function declaration for function to be integrated by ode45
function state_dot = odefun(~,state) 
    % Import global variables
    global g

    % Specify variable indexes within the initial condition array 
    x=state(1);
    y=state(2);
    z=state(3);
    xd=state(4);
    yd=state(5);
    zd=state(6);
    m=state(7);

    % Extract position and velocity vectors to be fed into acceleration
    % calculation
    v=[xd,yd,zd]';
    r=[x,y,z]';

    % Calculate acceleration at given time according using Optimal Guidance Law
    % paper derivations
    acceleration = tgolf(v, r, g);

    % Extract individual acceleration components
    xdd = acceleration(1);
    ydd = acceleration(2);
    zdd = acceleration(3);

    % Compute change in mass
    md = -norm(acceleration)*m/(9.81*237);

    state_dot=[xd,yd,zd,xdd, ydd, zdd, md]';
end 

%% Function declaration for event function to detect impact with ground
function [value,isterminal,direction] = Ground(~,state)
    %Event funtion to stop integration when vehicle reaches ground
    if state(3) > 1e-3
        value = 1; %Keep going
    else %If not
        value = 0; %Then stop
    end
    isterminal = 1; %Terminate integration when condtion met
    direction = 0; %Direction doesn't matter
end

