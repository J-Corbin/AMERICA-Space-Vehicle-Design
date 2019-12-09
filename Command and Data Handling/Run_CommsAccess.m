% This script is responsible for running the CommsAccess function over a
% series of ephemeris times, and tracking how the data transfer capacity
% from a lander to the orbiters changes over time. 
%
% Developed by William Parker, November 2019

clear all; close all; 
%set ephemeris times which exist within the windows of data transfer
%capability from the lander to the orbiters. change in increments of one
%martian sol (3600*24 seconds)
et1 = 536544069;
et2 = 599610000;
etOverall = et1:3600*24:et1+3600*24*687;

%Loop through each of the ephemeris times and record the Data_per_sol at
%each based on the orbiter configuration at that time. 
%     et = [epoch1:2592000:epoch2];
for i = 1:length(etOverall)
    [Data_per_sol(i)] = CommsAccess(etOverall(i)) ;
    i/687
%     plot((etOverall(i)-et1)/3600, Data_per_sol(i)*10^-9,'r');
% title('Data Capacity Per Sol over Martian Year');
% xlabel('Days');
% ylabel('Gb of Data per day')
% pause(0.02)
end
% 
% plot((etOverall-et1)/(3600*24), Data_per_sol/1e9,'r');
plot([1:length(etOverall)], Data_per_sol/1e9,'r');

% title('Data Transfer Capacity Per Sol over Martian Year');
xlabel('Sols');
ylabel('GBits of Data per Sol')
grid on 
