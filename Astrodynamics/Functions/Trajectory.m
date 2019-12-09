function [r, v, c3, v_inf_arrival, tof] = Trajectory(porkchop, launch_day, arrival_day, dark)
    %% Generate Pork Chop Plots
    % Generates porkchop plot if user requests
    if porkchop == 'y'
        Porkchop(launch_day, arrival_day)
    else
        fprintf('\n No porkchop plots will be displayed.\n')
        fprintf('------------------------------------------------------------------')
    end

    %% Constants Declaration
    % Gravitational parameters
    mu_sun = 1.32712e11;

    % Celestial body scale radii
    r_earth = .5e7;
    r_sun = 1.5e7; 
    r_mars = r_earth;

    % Colors for the celestial bodies
    yellow = [1 1 .1];
    blue = [0 .5 1];
    red = [.8 0 0];

    % Set intercept line to white if darkmode enabled
    if dark == 'y'
        intercept_color = 'w';
    else
        intercept_color = 'k';
    end

    %% State Vector Calculation
    % Calculate time (seconds) for launch, arrival, time of flight
    depart_time = cspice_str2et(launch_day);
    arrival_time = cspice_str2et(arrival_day);
    tof = arrival_time - depart_time;

    % Produce Earth and Mars state vectors for departure time
    [e_state_depart, ~] = cspice_spkezr('EARTH', depart_time, 'J2000', 'LT', 'SUN');
    [m_state_depart, ~] = cspice_spkezr('MARS', depart_time, 'J2000', 'LT', 'SUN');

    % Produce Earth and Mars state vectors for arrival time
    [e_state_arrival, ~] = cspice_spkezr('EARTH', arrival_time, 'J2000', 'LT', 'SUN');
    [m_state_arrival, ~] = cspice_spkezr('MARS', arrival_time, 'J2000', 'LT', 'SUN');

    %% Solving Lambert's Problem for Intercept Trajectory
    % Uncomment to use my lambert solver
    %[v1, v2] = lambert(e_state_depart(1:3), m_state_arrival(1:3), tof, mu_sun, 'pro');

    % Uncomment to use Dr. Christian's lambert solver
    [v1, v2] = glambert(mu_sun, [e_state_depart(1:3);e_state_depart(4:6)],...
        [m_state_arrival(1:3);m_state_arrival(4:6)], tof, 0);

    % Calculate delta-v at departure and arrival
    dv_depart = v1 - e_state_depart(4:6);
    dv_arrival = v2 - m_state_arrival(4:6);

    % Calculate departure C3 and v_inf at arrival
    c3 = dv_depart' * dv_depart;
    v_inf_arrival = norm(dv_arrival);

    % Output the state vectors of the final point of the interept trajectory
    r = m_state_arrival(1:3);
    v = v2;

    %% Calculating Max Distance Mars is from Earth During Trajectory
    % Specify time interval over which the transfer happens
    time_interval = depart_time:1000:arrival_time;
    
    % Initialize array to hold state vectors at each time
    e_state = zeros(6, length(time_interval));
    m_state = zeros(6, length(time_interval));
    
    % Initialize array to hold distances between both planets and satellite
    distance_planets = zeros(1, length(time_interval));
    distance_sat = zeros(1, length(time_interval));
    
    % Find the position of the satellite at each time over the trajectory
    [sat_pos, theta1, theta2] = sat_position(e_state_depart(1:3), v1, mu_sun,...
        m_state_arrival(1:3), v2, length(time_interval));
    % Find the state vectors for both planets at each point in time over
    % interval
    for i = 1:length(time_interval)
        % State vectors for Earth
       [e_state(:, i), ~] = cspice_spkezr('EARTH', time_interval(i),...
           'J2000', 'LT', 'SUN');
       % State vectors for Mars
       [m_state(:, i), ~] = cspice_spkezr('MARS', time_interval(i),...
           'J2000', 'LT', 'SUN');
       
       % Calculate distance between the two planets at each time
       distance_planets(i) = norm(e_state(1:3, i) - m_state(1:3, i));
       
       % Calculate the distance between Earth and the satellite at each
       % time
       distance_sat(i) = norm(e_state(1:3, i) - sat_pos(1:3, i));
    end
    
    % Display the max distance between the planets
    max_dist = max(distance_planets);
    fprintf('\n The max distance between the planets is %g km.\n', max_dist)
    fprintf('------------------------------------------------------------------')
    
    % Display the max distance between Earth and the spacecraft
    [max_dist, ~] = max(distance_sat);
    fprintf('\n The max distance between Earth and spacecraft is %g km.\n', max_dist)
    fprintf('------------------------------------------------------------------')
    
    %% Solar Conjunction Timing
    % Uncomment to perform solar conjunction analysis
    %{
    % Compute day the solar conjunction occurs
    [solar_conj_date, seconds] = solar_conjunction(launch_day);
    
    % Find position vectors to planets at conjunction
    [e_conj, ~] = cspice_spkezr('EARTH', seconds, 'J2000', 'LT', 'SUN');
    [m_conj, ~] = cspice_spkezr('MARS', seconds, 'J2000', 'LT', 'SUN');
    
    % Display the day of conjunction
    fprintf('\n The solar conjunction after departure occurs on: \n')
    fprintf("   " + solar_conj_date + '\n')
    fprintf('------------------------------------------------------------------')
   %}
    
    %% Plotting the Intercept Trajectory
    figure
    l_angle = theta1;
    u_angle = theta2;
    plot_orbit(m_state_arrival(1:3), v2, mu_sun, l_angle, u_angle, intercept_color, 'n');
    %line([e_state(1,index), sat_pos(1,index)], [e_state(2,index), sat_pos(2,index)], [e_state(3,index), sat_pos(3,index)], 'Color', 'k')

    %% Plotting Earth's Orbit
    l_angle = 0;
    u_angle = 2 * pi;
    plot_orbit(e_state_depart(1:3), e_state_depart(4:6), mu_sun, l_angle, u_angle, '--b', 'n')

    %% Plotting Mar's Orbit
    l_angle = 0;
    u_angle = 2 * pi;
    plot_orbit(m_state_depart(1:3), m_state_depart(4:6), mu_sun, l_angle, u_angle, '--r', 'n')

    %% Plotting Various Items
    % Plotting vector to positions and text, uncomment to show
    %plot3([0, e_state_depart(1)], [0, e_state_depart(2)], [0, e_state_depart(3)], '-.b')
    %plot3([0, m_state_depart(1)], [0, m_state_depart(2)], [0, m_state_depart(3)], '-.r')
    %plot3([0, e_state_arrival(1)], [0, e_state_arrival(2)], [0, e_state_arrival(3)], 'b')
    %plot3([0, m_state_arrival(1)], [0, m_state_arrival(2)], [0, m_state_arrival(3)], 'r')
    
    % Plotting coordinate system at Sun
    line([0, 1e8],[0, 0], [0, 0]), text(1.1e8, 0, 0, 'x')
    line([0, 0], [0, 1e8], [0, 0]), text(0, 1.1e8, 0, 'y')
    line([0, 0], [0, 0], [0 1e8]), text(0, 0, 1.1e8, 'z')
    
    % Plot the sun
    plot_planet(r_sun, yellow, [0, 0, 0])

    % Plot Earth and Mars at departure
    plot_planet(r_earth, blue, e_state_depart(1:3))
    plot_planet(r_mars, red, m_state_depart(1:3))

    % Plot Earth and Mars at arrival
    plot_planet(r_earth, blue, e_state_arrival(1:3))
    plot_planet(r_mars, red, m_state_arrival(1:3))

    % Uncomment to plot solar conjunction
    
    % Plot lines for solar conjunction
    %line([0, e_conj(1)], [0, e_conj(2)], [0, e_conj(3)], 'Color', 'black')
    %line([0, m_conj(1)], [0, m_conj(2)], [0, m_conj(3)], 'Color', 'black')
    %text(m_conj(1)*1.05, m_conj(2)*1.05, m_conj(3)*1.05, 'Solar Conjunction')
    
    
    % Random graph adjustments
    grid on
    axis equal
    xlabel('x [km]')
    ylabel('y [km]')
    zlabel('z [km]')
    view(3)
    
    %title('Orbits and Intercept Trajectory to Mars')
    legend('Intercept Trajectory', 'Earth Orbit', 'Mars Orbit', 'Location', 'best')

    % Print C3 and v_inf values to the command line
    fprintf('\n The C3 value for this launch is: %g km^2/s^2', c3)
    fprintf('\n The v_inf value at Mars is: %g km/s', v_inf_arrival)
    fprintf('\n')

    % Changing plot colors to DARKMODE if desired
    if dark == 'y'
        set(gca, 'color', [0 0 0])
        set(gca, 'xcolor', 'w')
        set(gca, 'ycolor', 'w')
        set(gca, 'zcolor', 'w')
        text(e_state_depart(1)*1.2, e_state_depart(2)*1.2, e_state_depart(3)*1.2, 'Departure', 'Color', 'W')
        text(m_state_depart(1)*1.2, m_state_depart(2)*1.2, m_state_depart(3)*1.2, 'Departure', 'Color', 'W')
        text(e_state_arrival(1)*1.2, e_state_arrival(2)*1.2, e_state_arrival(3)*1.2, 'Arrival', 'Color', 'W')
        text(m_state_arrival(1)*1.2, m_state_arrival(2)*1.2, m_state_arrival(3)*1.2, 'Arrival', 'Color', 'W')
    else
        text(e_state_depart(1)*1.2, e_state_depart(2)*1.2, e_state_depart(3)*1.2, 'Departure')
        %text(m_state_depart(1)*1.2, m_state_depart(2)*1.2, m_state_depart(3)*1.2, 'Departure')
        %text(e_state_arrival(1)*1.2, e_state_arrival(2)*1.2, e_state_arrival(3)*1.2, 'Arrival')
        text(m_state_arrival(1)*1.2, m_state_arrival(2)*1.2, m_state_arrival(3)*1.2, 'Arrival')
    end
    
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
end