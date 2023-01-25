clear

% Run learner solution.
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
Tmax = 10;
kappa = 0.1;
tau = 2;
nu  = 0.1;


% global a b c T N Vbar s0 r0 xi0 K Tmax kappa Gamma A 
N = size(b,2);
A = [1:N];
Gamma = [0,-1;1,0];

t = [0:kappa:Tmax];
r_set = [0:nu:T];

xout = xi0;

for i = 0:round(Tmax/tau)-1
    J = zeros(1, round(T/nu));
    for j = 0:round(T/nu)
        x0 = [xout(:, end); r_set(j+1); 0];
        [~,xi] = ode45(@(t,x) ode_sym(t,x, a, b, c, T, Vbar, K, Gamma, A), [0:kappa:Tmax], x0);
        J(j+1) = xi(end, 6);
    end
    r(i+1) = r_set(J == min(J));
    x0 = [xout(:, end); r(i+1); 0];
    [tiout,xiout]=ode45(@(t,x) ode_sym(t,x, a, b, c, T, Vbar, K, Gamma, A),[0:kappa:tau],x0);
    xout=[xout, xiout(2:end, 1:end-2)'];
end
s =2;
