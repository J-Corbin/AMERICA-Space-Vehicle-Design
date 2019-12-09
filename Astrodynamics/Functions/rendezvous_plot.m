function rendezvous_plot(soi_mars, periapsis, r_mars_atm, mu_mars, r_mars, v_inf)
    %Open a new figure
    figure

    %% Plot the SOI and the atmosphere
    %Calculate x, y, and z values for the atmosphere and SOI locations
    theta = 0:.01:2*pi;
    z = zeros(1, length(theta));
    %x_soi = soi_mars * cos(theta);
    %y_soi = soi_mars * sin(theta);
    x_atm = (r_mars_atm + r_mars) * cos(theta);
    y_atm = (r_mars_atm + r_mars) * sin(theta);
    
    %plot3(x_soi, y_soi, z, 'k')
    plot3(x_atm, y_atm, z, '--k')
    axis equal
    grid on
    hold on
    
    %% Plot Mars
    plot_planet(r_mars, [0.8 0 0], [0 0 0])
    
    %% Plot Trajectory
    %Calculate the delta and beta values for the approach trajectory
    delta = (periapsis + r_mars) * sqrt(1 + 2 * mu_mars / ...
        ((periapsis + r_mars) * v_inf^2));
    beta = asin(delta / soi_mars);
    
    %Giving state vectors a small z component to prevent division by 0
    r = zeros(3, 1) + 1e-8;
    v = zeros(3, 1) + 1e-8;
    r(1) = -soi_mars * cos(beta);
    r(2) = soi_mars * sin(beta);
    v(1) = v_inf;
    
    %Plot the approach trajectory to the atmosphere stopping at the
    %atmosphere
    [e, a, ~, ~, ~, ~] = orbit_elements(r, v, mu_mars);
    h = sqrt(a*mu_mars*(1 - e^2));
    angle = acos(h^2 / (mu_mars * (r_mars_atm + r_mars) * e) - 1 / e);
    plot_orbit(r, v, mu_mars, -pi/2, -angle, 'k', 'y')
    
    %Plot the portion of the trajectory that won't be travelled
    plot_orbit(r, v, mu_mars, -angle, pi/2, '--r', 'y')
    
    %% Plot Directional Arrows and Text
    
    line([0 5000], [0 0], [0 0], 'Color', 'b', 'LineWidth', 2)
    plot3(5000, 0, 0, 'b>')
    text(3750, -250, 0, 'Mars Velocity')
    line([0 0], [0 -5000], [0 0], 'Color', 'k', 'LineWidth', 2)
    plot3(0, -5000, 0, 'kv')
    text(400, -4750, 0, 'To Sun')
    
    
    %Figure adjustments
    legend('Edge of Atmosphere', 'Mars', 'Approach Trajectory (Path Taken)',...
        'Location of Periapsis', 'Departure Trajectory with no Atmosphere')
    xlabel('km')
    ylabel('km')
    zlabel('km')
    %title('Approach Trajectory to Mars')
    view(2)
    
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
end