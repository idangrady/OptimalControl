function xdot = ode_sym(t, x, a, b, c, T, V_bar, K, Gamma, A)
    % This function compares the desired position and current position and
    % applies feedback control to compensate for the error.

    % Compute intermediate values
    gamma = a + b * cos(2*pi/T*A'*x(5)) + c * sin(2*pi/T*A'*x(5));
    d_gamma = -b * 2*pi/T * (A' .* sin(2*pi/T*A'*x(5))) + c * 2*pi/T * (A' .* cos(2*pi/T*A'*x(5)));
    dd_gamma = -b * ((2*pi/T*A').^2 .* cos(2*pi/T*A'*x(5))) - c * ((2*pi/T*A').^2 .* sin(2*pi/T*A'*x(5)));

    % Compute desired values
    d_alpha = V_bar / sqrt(d_gamma'*d_gamma);
    xy_bar = gamma;
    dxy_bar = d_gamma * d_alpha;
    theta_bar = atan2(dxy_bar(2), dxy_bar(1));
    rho_bar = [cos(theta_bar); sin(theta_bar)];
    ddxy_bar = dd_gamma * (V_bar^2) / (d_gamma'*d_gamma) - d_gamma * (d_gamma'*dd_gamma) * V_bar^2 / (d_gamma'*d_gamma)^2;
    omega_bar = rho_bar' * Gamma' * ddxy_bar / V_bar;

    % Compute error values
    ex = x(1) - xy_bar(1);
    ey = x(2) - xy_bar(2);
    e_bar = [x(3), x(4); -x(4), x(3)] * [ex; ey];
    ex_bar = e_bar(1);
    ey_bar = e_bar(2);
    e_rho = 1 - x(3:4)' * rho_bar;

    % Compute h value
    if x(3:4)' * rho_bar >= 0
        h = -x(3:4)' * Gamma' * rho_bar;
    elseif x(3:4)' * rho_bar < 0 && x(3:4)' * Gamma' * rho_bar < 0
        h = 1;
    else
        h = -1;
    end

    % Compute control inputs
    V = V_bar * (1 - e_rho) - K(1) * ex_bar;
    omega = omega_bar - K(2) * V_bar * ey_bar - K(3) * h;

    % Compute state derivatives
    dx = V * x(3);
    dy = V * x(4);
    d_rho = [0,-omega; omega, 0] * x(3:4);
xdot = [dx; dy; d_rho; d_alpha];
end