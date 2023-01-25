function x_dot = ode_sym(t, x, a, b, c, T, V_des, K, Gamma, A)
    % Define intermediate variables
    gamma = a + b * cos(2 * pi / T * A' * x(5)) + c * sin(2 * pi / T * A' * x(5));
    d_gamma = - b * 2 * pi / T * (A' .* sin(2 * pi / T * A' * x(5))) + c * 2 * pi / T * (A' .* cos(2 * pi / T * A' * x(5)));
    dd_gamma = - b * ((2 * pi / T * A').^2 .* cos(2 * pi / T * A' * x(5))) - c * ((2 * pi / T * A').^2 .* sin(2 * pi / T * A' * x(5)));
    d_alpha = V_des / norm(d_gamma);
    xy_des = gamma;
    dx_des = d_gamma * d_alpha;
    theta_des = atan2(dx_des(2), dx_des(1));
    rho_des = [cos(theta_des); sin(theta_des)];
    ddx_des = dd_gamma * (V_des^2) / (d_gamma' * d_gamma) - d_gamma * dot(d_gamma, dd_gamma) * V_des^2 / (d_gamma' * d_gamma)^2;
    omega_des = rho_des' * Gamma' * ddx_des / V_des;


    % Define error variables
    e_x = x(1) - xy_des(1); 
    e_y = x(2) - xy_des(2);
    e_bar = [x(3), x(4); -x(4), x(3)] * [e_x; e_y];
    e_x_bar = e_bar(1); e_y_bar = e_bar(2);
    e_rho = 1 - x(3:4)' * rho_des;

    % Determine h value
    switch sign(x(3:4)' * rho_des)
        case 1
            h = -x(3:4)' * Gamma' * rho_des;
        case 0
            h = 1;
        case -1
            h = -1;
    end     

    % Define state derivatives
    V = V_des * (1 - e_rho) - K(1) * e_x_bar;
    omega = omega_des - K(2) * V_des * e_y_bar - K(3) * h;
    dx = V * x(3); dy = V * x(4);
    drho = [0, -omega; omega, 0] * x(3:4);
    dJ = e_x * e_x + e_y * e_y + e_rho;
    x_dot = [dx; dy;drho; d_alpha;dJ];
end
