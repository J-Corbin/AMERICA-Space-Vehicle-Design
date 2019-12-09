function Porkchop(launch_day, arrival_day)
    %%%%%%%%%%%%%%%%
    %
    % Porkchop.m
    % This code generates porkchop plot contours for interplanetary transfers.
    % Default is Earth-Mars transfer for 2005 opportunity.
    %
    % Requies the following toolboxes:
    %   >> NASA JPL SPICE Tookit for MATLAB (MICE):
    %           http://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html
    %   >> Lambert's Problem Toolbox (from MATLAB File Exchange, by David Eagle)
    %           http://www.mathworks.com/matlabcentral/fileexchange/39530-lambert-s-problem
    %
    % Primary Developer Contact Information:
    % Dr. John A. Christian
    % Assistant Professor, Dept. of Mechanical, Aerospace, and Nuclear Engineering
    % Rensselaer Polytechnic Institute
    % chrisj9@rpi.edu
    %
    % Development History
    % Date       Developer       Comments
    % --------   -------------   ----------------------
    % 03/09/15   J. Christian    Initial implementation
    % 09/12/19   J. Christian    Updated function calls
    % 10/3/19    J. Corbin       Restructured code to be a function
    %
    % (c) Dr. John Christian, 2019
    %
    %%%%%%%%%%%%%%%%
    day2sec = 24*3600;
    sec2day = 1/day2sec;

    cspice_furnsh('naif0012.tls');
    cspice_furnsh('de421.bsp');
    cspice_furnsh('pck00010.tpc');

    %Bounds for launch/arrival plots
    t_depart='Mar 18, 2024 12:00:00.000000';
    t_arrive='Dec 18, 2024 12:00:00.000000';
    
    %There are 2099 days between Jan 1, 2020 and Sep 30, 2025, our NLT date.
    window_depart_days = 500;
    window_arrive_days = 800;

    %Find departure time (et is in units of seconds)
    et_depart=cspice_str2et(t_depart);

    %Find arrival time (et is in units of seconds)
    et_arrive=cspice_str2et(t_arrive);

    %Actual launch/arrival dates
    %launch_day = 'Sep 10, 2022 12:00:00.000000';
    depart = (cspice_str2et(launch_day) - et_depart) * sec2day;
    %arrival_day = 'Mar 28, 2023 12:00:00.000000';
    arrival = (cspice_str2et(arrival_day) - et_arrive) * sec2day;
    %time_orbit = (cspice_str2et(arrival_day) - cspice_str2et(launch_day)) * sec2day;

    SweepStepDepart = 300;
    SweepStepArrive = 300;
    et_depart_sweep = linspace(et_depart,et_depart+window_depart_days*day2sec, SweepStepDepart);
    et_arrive_sweep = linspace(et_arrive,et_arrive+window_arrive_days*day2sec, SweepStepArrive);

    GM_SUN = 132712441933.0;

    C3_depart   = zeros(SweepStepDepart,SweepStepArrive);
    vinf_arrive = zeros(SweepStepDepart,SweepStepArrive);

    for ii = 1:SweepStepDepart

        %Get Earth position and velocity
        et_ii = et_depart_sweep(ii);
        [state, ~] = cspice_spkezr('EARTH', et_ii, 'J2000', 'LT', 'SUN');
        r_earth = state(1:3);
        v_earth = state(4:6);

        for jj= 1:SweepStepArrive

            %Get Mars position and velocity
            et_jj = et_arrive_sweep(jj);
            [state, ~] = cspice_spkezr('Mars', et_jj, 'J2000', 'LT', 'SUN');
            r_mars = state(1:3);
            v_mars = state(4:6);

            %Compute time of flight (in seconds)
            tf = et_jj - et_ii;

            %Limit search to reasonable TOFs
            if tf*sec2day >= 100 && tf*sec2day <= 700

                %Solve Lambert's problem
                [vi, vf] = glambert(GM_SUN, [r_earth;v_earth], [r_mars;v_mars], tf, 0);
                dv_depart = vi-v_earth;
                dv_arrive = vf-v_mars;

                %Compute departure C3 and arrival Vinf
                C3_depart(ii,jj) = dv_depart'*dv_depart;
                vinf_arrive(ii,jj) = norm(dv_arrive);
            else
                C3_depart(ii,jj) = inf;
                vinf_arrive(ii,jj) = inf;
            end

        end
    end

    %% Departure Figure
    levels = 0:.05:30;
    figure
    [~,h] = contour( (et_depart_sweep-et_depart)*sec2day,(et_arrive_sweep-et_arrive)...
        *sec2day,C3_depart','ShowText','off');
    axis equal
    grid on
    set(h,'LevelList',levels)
    xlabel('Days after March 18, 2024 (Departure)')
    ylabel('Days after December 18, 2024 (Arrival)')
    %title('Contours of Departure $C3 = V_{inf}^2$ ($\frac{km^2}{s^2}$)', 'Interpreter', 'latex')
    %xtickangle(-45)
    %xticklabels({'Mar 18, 2022 - 0', 'Jun 26, 2022 - 100', 'Oct 4, 2022 - 200',...
    %    'Jan 1, 2023 - 300', 'Apr 4, 2023 - 400', 'Jul 31, 2023 - 500'}) 
    %ytickangle(45)
    %yticklabels({'Dec 18, 2022 - 0', 'Mar 28, 2023 - 100', 'Jul 6, 2023 - 200',...
    %    'Oct 14, 2023 - 300', 'Jan 22, 2024 - 400', 'May 1, 2024 - 500', ...
    %    'Aug 9, 2024 - 600', 'Nov 17, 2024 - 700', 'Feb 25, 2025 - 800'})

    %Plotting line for launch
    line([depart, depart], [0, 800], 'Color', 'red', 'LineStyle', '--')
    line([0, 500], [arrival, arrival], 'Color', 'red', 'LineStyle', '--')
    hold on
    plot(depart, arrival, 'r*')
    b = colorbar;
    b.Label.String = 'km^2/s^2';
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    
    %% Arrival Figure
    figure
    levels = 0:.01:15;
    [~,h] = contour( (et_depart_sweep-et_depart)*sec2day,(et_arrive_sweep-et_arrive)...
        *sec2day,vinf_arrive','ShowText','off');
    axis equal
    grid on
    set(h,'LevelList',levels)
    xlabel('Days after March 18, 2024 (Departure)')
    ylabel('Days after December 18, 2024 (Arrival)')
    %title('Contours of Arrival $V_{inf}$ ($\frac{km}{s}$)', 'Interpreter', 'latex')
    %xtickangle(-45)
    %xticklabels({'Mar 18, 2022 - 0', 'Jun 26, 2022 - 100', 'Oct 4, 2022 - 200',...
    %    'Jan 1, 2023 - 300', 'Apr 4, 2023 - 400', 'Jul 31, 2023 - 500'}) 
    %ytickangle(45)
    %yticklabels({'Dec 18, 2022 - 0', 'Mar 28, 2023 - 100', 'Jul 6, 2023 - 200',...
    %    'Oct 14, 2023 - 300', 'Jan 22, 2024 - 400', 'May 1, 2024 - 500',...
    %    'Aug 9, 2024 - 600', 'Nov 17, 2024 - 700', 'Feb 25, 2025 - 800'})

    %Plotting line for launch
    line([depart, depart], [0, 800], 'Color', 'red', 'LineStyle', '--')
    line([0, 500], [arrival, arrival], 'Color', 'red', 'LineStyle', '--')
    hold on
    plot(depart, arrival, 'r*')
    a = colorbar;
    a.Label.String = 'km/s';
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
end