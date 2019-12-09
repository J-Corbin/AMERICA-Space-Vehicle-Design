function plot_planet(radius, color, location)
    %% Draw Planet
    [xx, yy, zz] = sphere(100); %Creates sphere coordinates
    surface(radius*xx + location(1), radius*yy + location(2),...
        radius*zz + location(3), 'FaceColor', color, 'EdgeColor', 'none',...
        'FaceAlpha', 0.5);
        %Draw the sphere as a surface
end