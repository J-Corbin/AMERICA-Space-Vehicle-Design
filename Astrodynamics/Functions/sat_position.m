function [coords, theta1, theta2] = sat_position(r, v, u, r2, v2, samples)
    deg_rad = pi / 180;

    [mag_e, ~, ~, Omega, inc, omega] = orbit_elements(r, v, u);

    %Magnitude of momentum for unknown object
    h = norm(cross(r, v));
    
    % Finding lower angle limit
    [~, ~, theta1, ~, ~, ~] = orbit_elements(r, v, u);
    theta1 = deg_rad * theta1;
    
    % Finding upper angle limit
    [~, ~, theta2, ~, ~, ~] = orbit_elements(r2, v2, u);
    theta2 = deg_rad * theta2;
    
    % Determine step size
    step = (theta2 - theta1) / (samples - 1);
    
    %Converting estimated orbital parameters to radians
    W = Omega * deg_rad;
    w = omega * deg_rad;
    inclination = inc * deg_rad;

    %Computing direction cosine matrix from perifocal to N for unknown object
    C_N_PF_obj = [-sin(W)*cos(inclination)*sin(w)+cos(W)*cos(w), -sin(W)*cos(inclination)*cos(w)-cos(W)*sin(w), sin(W)*sin(inclination);...
    cos(W)*cos(inclination)*sin(w)+sin(W)*cos(w), cos(W)*cos(inclination)*cos(w) - sin(W)*sin(w), sin(inclination)*cos(w);...
    sin(inclination)*sin(w), sin(inclination)*cos(w), cos(inclination)];

    %For loop to calculate X,Y,Z positions of the estimated orbit
    counter = 1;
    for angle = theta1:step:theta2
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

    coords = coordinate;
end