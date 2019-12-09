function plot_orbit(r, v, u, l_limit, u_limit, color, draw)
    deg_rad = pi / 180;

    [mag_e, ~, ~, Omega, inc, omega] = orbit_elements(r, v, u);

    %Converting estimated orbital parameters to radians
    W = Omega * deg_rad;
    w = omega * deg_rad;
    inclination = inc * deg_rad;

    %Computing direction cosine matrix from perifocal to N for unknown object
    C_N_PF_obj = [-sin(W)*cos(inclination)*sin(w)+cos(W)*cos(w), -sin(W)*cos(inclination)*cos(w)-cos(W)*sin(w), sin(W)*sin(inclination);...
    cos(W)*cos(inclination)*sin(w)+sin(W)*cos(w), cos(W)*cos(inclination)*cos(w) - sin(W)*sin(w), sin(inclination)*cos(w);...
    sin(inclination)*sin(w), sin(inclination)*cos(w), cos(inclination)];

    %Magnitude of momentum for unknown object
    h = norm(cross(r, v));

    %For loop to calculate X,Y,Z positions of the estimated orbit
    counter = 1;
    for angle = l_limit:.001:u_limit
    radius_polar = [h^2/(u * (1 + mag_e * cos(angle)));0;0];
    C_PF_polar_obj = [cos(angle), sin(angle), 0;...
    -sin(angle), cos(angle), 0;...
    0, 0, 1].';
    radius_pf = (C_PF_polar_obj * radius_polar); 
    radius = (C_N_PF_obj * radius_pf).';
    for z = 1:3
       coordinate(z, counter) = radius(z);
    end
    counter = counter + 1;
    end
    per = 0;

    radius_polar = [h^2/(u * (1 + mag_e * cos(per)));0;0];
    C_PF_polar_obj = [cos(per), sin(per), 0;...
    -sin(per), cos(per), 0;...
    0, 0, 1].';
    radius_pf = (C_PF_polar_obj * radius_polar); 
    radius = (C_N_PF_obj * radius_pf).';
    periapsis = radius;

    %Plot the previously determined orbit
    plot3(coordinate(1,:), coordinate(2,:), coordinate(3,:), color)
    hold on

    %Draw line to periapsis and a dot at periapsis
    if draw == 'y'
        %line([0 periapsis(1)], [0 periapsis(2)], [0 periapsis(3)],...
            %'Color', 'k', 'LineStyle', '--')
        text(periapsis(1)*1.1, periapsis(2)*1.1, periapsis(3)*1.1, 'periapsis')
        %plot_planet(100, 'k', periapsis)
        plot3(periapsis(1), periapsis(2), periapsis(3), 'k*')
    end
end