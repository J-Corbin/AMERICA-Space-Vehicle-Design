function [date, et] = solar_conjunction(start_time)
    % Convert times to seconds
    time_count = cspice_str2et(start_time);
    diff = time_count;
    year_to_second = 365.25 * 24 * 60 * 60;
    years = 3;

    % Constantly iterate until solar conjunction is found
    while true
        % Determine position of planets
        [e_state, ~] = cspice_spkezr('EARTH', time_count, 'J2000', 'LT', 'SUN');
        [m_state, ~] = cspice_spkezr('MARS', time_count, 'J2000', 'LT', 'SUN');

        % Create unit vectors of planet positions
        e_unit_vec = e_state(1:3) ./ norm(e_state(1:3));
        m_unit_vec = m_state(1:3) ./ norm(m_state(1:3));
        e_unit = round(e_unit_vec, 2);
        m_unit = round(m_unit_vec, 2);
        
        % Add the two unit vectors
        sum = e_unit + m_unit;

        % Check if the magnitude of the unit vector is below a threshold
        if norm(sum) < 0.05
            break
        end
        
        % End loop if time is exceeded
        if time_count - diff >= years * year_to_second
            break
        end
        
        % Iterate to next timestamp
        time_count = time_count + 10;
    end
    
    % Format how time should be presented
    sample = 'Thu Oct 1 11:11:11 PDT 1111';
    [pictur, ok, errmsg] = cspice_tpictr(sample);
    
    % Save ephemeris time value
    et = time_count;
    
    % Save date of conjunction
    date = cspice_timout(time_count, pictur);
end