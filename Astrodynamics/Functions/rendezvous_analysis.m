function rendezvous_analysis(soi, altitude, atmosphere, mu, radius, v_inf)
    %% Set range of periapsis values and calculate resultant geometry
    r_periapsis = 1:(radius + altitude); %Array of various periapsis radii
    delta = r_periapsis .* sqrt(1 + 2 * mu ./...
        (r_periapsis .* v_inf^2)); %Aiming radii 
    beta = asin(delta ./ soi); %Angle from Mars horizontal to delta at the
    %SOI

    % Initialize arrays to hold various parameters
    r = zeros(3, length(r_periapsis)) + 1e-8;
    v = zeros(3, length(r_periapsis)) + 1e-8;
    h = zeros(1, length(r_periapsis)) + 1e-8;
    e = zeros(1, length(r_periapsis)) + 1e-8;
    
    %Filling the previously initialized arrays for position and velocity
    r(1,:) = -soi * cos(beta);
    r(2,:) = soi * sin(beta);
    v(1,:) = v_inf;

    %Calculating momentum and eccentricity
    for i = 1:length(r_periapsis)
        h(i) = norm(cross(r(:,i), v(:,i)));
        e(i) = norm(cross(v(:,i), cross(r(:,i), v(:,i))) / mu - ...
            r(:,i) / norm(r(:,i)));
    end

    %Calculating anomally where orbit intersects atmosphere and flight path
    %angle at that point
    anomally = acos((h.^2 / (mu * (atmosphere + radius)) - 1) ./ e);
    flight_path_angle = atan((e .* sin(anomally)) ./ (1 + e .*...
        cos(anomally)));

    %Plot computed data for flight path angle as function of periapsis
    figure
    plot(r_periapsis, flight_path_angle * 180 / pi, 'k')
    grid on
    title('Flight Path Angle')
    xlabel('periapsis (km)')
    ylabel('Flight Path Angle (Deg)')
    xlim([0 radius + altitude])
    line([radius radius], [0 90], 'Color', 'red', 'LineStyle', '--')
    legend('Flight Path Angle', 'Mars Surface','Location', 'southwest')

    %Calculate the atmospheric entry velocity
    energy = v_inf^2 / 2;
    v_atm = sqrt(2 * (energy + mu / (radius + altitude)));
    fprintf('The velocity at atmosphere is: %g km/s\n', v_atm)

    %Uncomment to print data to a .txt file
    %{
    fileID = fopen('fpa.txt', 'w');
    fprintf(fileID, '%s %s\n', 'periapsis (km)', 'Flight Path Angle (deg)');
    fprintf(fileID,'%6.2f %12.8f\n',[r_periapsis; flight_path_angle * 180 / pi]);
    fclose(fileID);
    %}
end