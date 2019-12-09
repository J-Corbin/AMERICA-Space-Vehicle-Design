%% Plotting tool for Power Profiles
%  Power profile of surface operations during cycle of 7 sols

% pulling data from spreadsheet
time  = transpose(xlsread('schedule.xlsx','Sheet1','a2:a673'));
power = transpose(xlsread('schedule.xlsx','Sheet1','b2:b673'));

% plotting
plot(time,power)
xlabel('Time (h)')
ylabel('Power (W)')
axis([0 167.75 0 350])
grid on
hold off