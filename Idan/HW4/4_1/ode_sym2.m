function xdot = ode_sym2(t, x, a, b, c, T, V_bar, K, Gamma, A)
    % Compare desired and current position and apply feedback control
    
    % Intermediate values
    gamma = a + b*cos(2*pi/T*A'*x(5)) + c*sin(2*pi/T*A'*x(5)); 
    dgamma = -b*2*pi/T*(A'.*sin(2*pi/T*A'*x(5))) + c*2*pi/T*(A'.*cos(2*pi/T*A'*x(5)));  
    ddgamma = -b*((2*pi/T*A').^2.*cos(2*pi/T*A'*x(5))) - c*((2*pi/T*A').^2.*sin(2*pi/T*A'*x(5)));
    
    % Desired values
    dalpha = V_bar/sqrt(dgamma'*dgamma);
    xy_bar = gamma;
    dxy_bar = dgamma*dalpha;
    theta_bar = atan2(dxy_bar(2), dxy_bar(1));
    rhobar = [cos(theta_bar);sin(theta_bar)];
    ddxy_bar = ddgamma*(V_bar^2)/(dgamma'*dgamma) - dgamma*dgamma'*ddgamma*V_bar^2/(dgamma'*dgamma)^2;
    omegabar = rhobar'*Gamma'*ddxy_bar/V_bar;

    % Error values
    ex = x(1) - xy_bar(1);
    ey = x(2) - xy_bar(2);
    ebar = [x(3), x(4); -x(4), x(3)]*[ex;ey];
    exbar = ebar(1);
    eybar = ebar(2);
    erho = 1 - x(3:4)'*rhobar;
    
    % Control inputs
    V = V_bar*(1 - erho) - K(1)*exbar;
    if x(3:4)'*rhobar >= 0
        h = -x(3:4)'*Gamma'*rhobar;
    elseif x(3:4)'*rhobar < 0 && x(3:4)'*Gamma'*rhobar < 0
h = 1;
else %if x(3:4)'*rhobar < 0 && x(3:4)'*Gamma'*rhobar>0
h = -1;
end
omega = omegabar - K(2)*eybar - K(3)*h;
% State derivatives
dx = V*x(3);
dy = V*x(4);
drho = [0, -omega; omega, 0]*x(3:4);
xdot = [dx; dy; drho; dalpha];
end