clear, clc, close all

%% NASA Data
% Load the NASA planetary data modules
cspice_furnsh('naif0012.tls');
cspice_furnsh('de421.bsp');
cspice_furnsh('pck00010.tpc');
    
%% Variable Declaration
% Gravitational Parameters
mu_mars = 42828;
mu_sun = 1.32712e11;

% Shape of Mars
r_mars = 3389.5; %Mars' radius (km)
r_mars_atm = 125; %Height of Mars' atmosphere (km)
soi_mars = 577000; %Radius of Mars' sphere of influence (km)

% Unit conversion
deg_rad = pi / 180; %converting degree to rad
day2sec = 24*3600;
sec2day = 1/day2sec;

%Earlier dates: launch: Sep 10, 2022
%               arrive: Apr 16, 2023
%Later dates: Launch: Oct 25 10, 2024
%               arrive: Jun 11, 2025

launch_day = 'Oct 25, 2024 12:00:00.000000';
arrival_day = 'Jun 11, 2025 12:00:00.000000';

%% Transfer Trajectory
% Plot the orbits of Earth and Mars, along with the transfer trajectory
[r, v, c3, v_inf, tof] = Trajectory('y', launch_day, arrival_day, 'n');
%close all %Comment to show figure

% Calculating time of flight for the transfer
tof = tof * sec2day;
fprintf('------------------------------------------------------------------')
fprintf('\n The transfer to Mars will take approximately %g days.\n', tof)

% Calculate the orbital elements of the transfer trajectory
[mag_e, a, ~, Omega, i, omega] = orbit_elements(r, v, mu_sun);
fprintf('------------------------------------------------------------------')
fprintf('\n The orbital elements of the transfer trajectory are:')
fprintf('\n Eccentricity =                                       %g', mag_e)
fprintf('\n Semi-Major Axis [km] =                               %g', a)
fprintf('\n Right of Ascension of the Ascending Node [degrees] = %g', Omega)
fprintf('\n Inclination [degrees] =                              %g', i)
fprintf('\n Argument of periapsis [degrees] =                    %g', omega)
fprintf('\n')
fprintf('------------------------------------------------------------------')
fprintf('\n')

%% Mars Rendezvous - Analysis
%Set the upper limit for periapsis analysis
max_periapsis_altitude = .85 * r_mars_atm;

%Throws error if max periapsis is outside the atmosphere
if max_periapsis_altitude > r_mars_atm
    clear, clc, close all
    error('The desired max periapsis altitude for analysis must be inside the atmosphere.')
end

%Determine if user wants to run flight path analysis
query = input('Do you want to run the flight path angle analysis? (y/n)\n', 's');
%fprintf('No flight path analysis. \n')
%query = 'y';
if query == 'y'
    %Analyze flight path angles and v_atm values for a range of periapsiss
    rendezvous_analysis(soi_mars, max_periapsis_altitude, r_mars_atm, mu_mars, r_mars, v_inf);
end
fprintf('------------------------------------------------------------------\n')

%% Mars Rendezvous - Plotting
%Constants determined from rendezvous analysis
periapsis = 3366 - 3389.5; %height of periapsis above Mars' surface (km)

%Plot the rendezvous trajectory
rendezvous_plot(soi_mars, periapsis, r_mars_atm, mu_mars, r_mars, v_inf)