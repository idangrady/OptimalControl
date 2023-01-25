clear all; close all; clc;
a    = [1;1];
b    = [0 0 0;1 0 0];
c    = -[1 0 0;0 0 0];
T    = 10;
N    = size(b,2);
Vbar = 1;
s0   = 0;
r0   = 0;
xi0  = [1;1;1;0];
K    = [1 1 1];
Tmax = 30;
kappa = 0.1;

% global a b c T N Vbar s0 r0 xi0 K Tmax kappa Gamma A 
N    = size(b,2);
A=[1:N];
Gamma=[0,-1;1,0];

t=[0:kappa:Tmax];
alpha0=r0; %Initial Condition for alpha(given in the document)
x0=[xi0;alpha0]; %Initial Condition
[ti1,xi1]=ode45(@(t,x) ode_sym(t,x, a, b, c, T, Vbar, K, Gamma, A),[0:kappa:Tmax],x0);
[ti,xi]=ode45(@(t,x) ode_sym2(t,x, a, b, c, T, Vbar, K, Gamma, A),[0:kappa:Tmax],x0);

xi = interp1(ti,xi,t); %Adjust the timestep between data points to kappa
x=xi(:, 1)'; y=xi(:, 2)'; rho=xi(:, 3:4)';