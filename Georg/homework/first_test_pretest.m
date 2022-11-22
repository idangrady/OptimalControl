omegai = [0 50 0];
tau    = 0.2;
Mv      = 30;
L2      = 10;
Na      = 5;
deltav  = 0.05;
si     = [0.3 0.4 floor(0.645/(deltav*tau))*deltav*tau];

[c] = timeoptimalpathspeed(omegai,si,tau,deltav,L2,Na,Mv);
% Run reference solution.
