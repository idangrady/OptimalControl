%Given Parameters
clear all; close all; clc;
a    = [0;1];
b    = [0 0 0;1 0 0];
c    = -[1 0 0;0 0 0];
T    = 10;
N    = size(b,2);
Vbar = 2;
s0   = 0;
r0   = 0;
xi0  = [1;1;1;0];
K    = [2 0.5 0.5];
Tmax = 30;
kappa = 0.1;
tau = 2;
nu  = 0.1;
%% Preparation for running ODE45
% global a b c T N Vbar s0 r0 xi0 K Tmax kappa Gamma A 
N    = size(b,2);
A=[1:N];
Gamma=[0,-1;1,0];
r_sym = sym('r_sym');
t=[0:kappa:Tmax];
r_set=[0:nu:T];
t_tau=[0:tau:Tmax];
%%
% Compute the optimal set of values r which minimize the cost function with the given horizon H=T
% using Dynamic Programming in the inner loop. 
% Once r is found, use r(1) to update the output values in the outer loop. 
tout=0;
xout=xi0;
J=0;
for i=0:round(Tmax/tau)-1
    for j=0:round(T/nu)
        x0=[xout(:, end);r_set(j+1);0]; %Initial Condition
        [ti,xi]=ode45(@(t,x) ode_sym(t,x, a, b, c, T, Vbar, K, Gamma, A),[0:kappa:Tmax],x0);
        J(j+1)=xi(size(xi, 1), 6);
    end
    r(i+1)=r_set(J==min(J));
    x0=[xout(:, end);r(i+1);0]; %Initial Condition
    [tiout,xiout]=ode45(@(t,x) ode_sym(t,x, a, b, c, T, Vbar, K, Gamma, A),[0:kappa:tau],x0);
    xout=[xout, xiout(2:end, 1:end-2)'];
%     tout=[tout, tiout(2:end)'];
end


% alpha0=1; %Initial Condition for alpha(given in the document)
% x0=[xi0;alpha0]; %Initial Condition
% [ti,xi]=ode45(@(t,x) ode_sym(t,x, a, b, c, T, Vbar, K, Gamma, A),[0,Tmax],x0);
% xout = interp1(tout',xout',t); %Adjust the timestep between data points to kappa
% set_des = interp1(t_tau,[x_des, y_des, rho_des],t); %Adjust the timestep between data points to kappa
x=xout(1, :); y=xout(2, :); rho=xout(3:4, :);
% xout=xout';
% x_des=set_des(:, 1)'; y_des=set_des(:, 2)'; rho_des=set_des(:, 3:4)';
%% Visualization
figure()
plot(x, y)

%%
function xdot=ode_sym(t,x, a, b, c, T, Vbar, K, Gamma, A) 
%This ode function compares the desired position and current position and
%Feedback controller is implemented to compensate for the erorr. All the
%equations are in the document variables with 'bar' are the desired values.
%     home, t %This is here to monitor the progress of the ode function. Remove or Comment out this line when running in Matlab Grader.
    
    gamma=a+b*cos(2*pi/T*A'*x(5))+c*sin(2*pi/T*A'*x(5)); 
    dgamma=-b*2*pi/T*(A'.*sin(2*pi/T*A'*x(5)))+c*2*pi/T*(A'.*cos(2*pi/T*A'*x(5)));  
    ddgamma=-b*((2*pi/T*A').^2.*cos(2*pi/T*A'*x(5)))-c*((2*pi/T*A').^2.*sin(2*pi/T*A'*x(5)));
    dalpha=Vbar/sqrt(transpose(dgamma)*dgamma);
    xybar=gamma;
    dxybar=dgamma*dalpha;
    thetabar=atan2(dxybar(2), dxybar(1));
    rhobar=[cos(thetabar);sin(thetabar)];
    ddxybar=ddgamma*(Vbar^2)/(transpose(dgamma)*dgamma)-dgamma*dot(dgamma, ddgamma)*Vbar^2/(transpose(dgamma)*dgamma).^2;
    omegabar=transpose(rhobar)*transpose(Gamma)*ddxybar/Vbar;

    ex=x(1)-xybar(1); 
    ey=x(2)-xybar(2);
    
    ebar=[x(3), x(4); -x(4), x(3)]*[ex;ey];
    
    exbar=ebar(1); eybar=ebar(2);
    
    erho=1-x(3:4)'*rhobar;
    
    if x(3:4)'*rhobar >=0
        h=-x(3:4)'*Gamma'*rhobar;
    elseif x(3:4)'*rhobar < 0 && x(3:4)'*Gamma'*rhobar<0
        h=1;
    else %if x(3:4)'*rhobar < 0 && x(3:4)'*Gamma'*rhobar>0
        h=-1;
    end     
    V=Vbar*(1-erho)-K(1)*exbar;
    omega=omegabar-K(2)*Vbar*eybar-K(3)*h;
    dx=V*x(3); dy=V*x(4);
    drho=[0, -omega;omega, 0]*x(3:4);
    dJ=ex*ex+ey*ey+erho;
    xdot=[dx;dy;drho;dalpha;dJ];
end