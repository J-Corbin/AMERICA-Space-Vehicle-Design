function [posN,t] = MyOrbitWithAnomaly(t,mu,a,e,h,incl,RAAN,ArgPer,anomaly)
%% Function Description
% This function is responsible for taking orbit parameters (including actual 
% anomaly) as inputs, and outputting a position vector at each time interval t 
% from the center of the attracting body to the object of interest based on the
% orbit parameters. 
% 
% INPUTS: 
% t = nx1 array of times since initial actual anomaly
% mu = gravitational parameter of the attracting body
% a = semimajor axis 
% e = eccentricity of orbit
% h = angular momentum
% incl = inclination (degrees)
% RAAN = right ascention of the ascending node (degrees)
% ArgPer = argument of periapsis (degrees)
% anomaly = actual anomally at the first ephemeris time of measurement
% OUUTPUTS:
% PosN = an nx3 matrix of position vectors with 3 components for each time t
% t = nx1 array of times since initial actual anomaly
% 
% Developed by William Parker, November 2019

%% Initialize variables
syms E theta
solE = zeros(length(t),1);
theta = zeros(length(t),1);
r = zeros(length(t),1);
PosOrb = zeros(3,length(t));
posPF = zeros(3,length(t));

%Create direction cosine matrix for converting a newtonian frame to the
%perifocal
C_N2PF = [-sin(RAAN)*cos(incl)*sin(ArgPer)+cos(RAAN)*cos(ArgPer), cos(RAAN)*cos(incl)*sin(ArgPer)+sin(RAAN)*cos(ArgPer), sin(incl)*sin(ArgPer);...
            -sin(RAAN)*cos(incl)*cos(ArgPer)-cos(RAAN)*sin(ArgPer), cos(RAAN)*cos(incl)*cos(ArgPer)-sin(RAAN)*sin(ArgPer), sin(incl)*cos(ArgPer);...
            sin(RAAN)*sin(incl), -cos(RAAN)*sin(incl), cos(incl)];

%for each time t, use Kepeler's equation to solve for orbital
%position , then convert to perifocal, then to Newtonian
for i = 1:length(t)
% Kep(i) = sqrt(mu/a)*(t(i)) == e*sin(E)-E;
solE(i) = vpasolve(sqrt(mu/a^3)*(t(i)) == E-e*sin(E),E);
theta(i) = 2*atan(sqrt((1+e)/(1-e))*tan(solE(i)/2))+anomaly;
% solTheta(i) = vpasolve(tan(theta/2)==sqrt((1+e)/(1-e))*tan(solE(i)/2),theta);
r(i) = h^2/mu*(1+e*cos(theta(i))).^-1;
PosOrb(:,i) = [r(i);0;0];
C_orb2PF = [cos(theta(i)) sin(theta(i)) 0;-sin(theta(i)) cos(theta(i)) 0;0 0 1];
posPF(:,i) = ((C_orb2PF)'*PosOrb(:,i));
posN(:,i) = C_N2PF'*posPF(:,i);
end

end