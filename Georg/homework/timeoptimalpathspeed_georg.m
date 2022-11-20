function [c] = timeoptimalpathspeed(omegai,si,tau,deltav,L2,Na,Mv)
%parameters:
%omegai     vector,   given omegavector
%si...      vector,   given distance vector
%tau...     constant, stepsize for the time
%deltav...  constant, stepsize for speed
%L2...      constant, limit for velocity change
%Na...      constant
%Mv...      constant


%1st calculate limits
L1 = deltav*Na/tau;
%L2 is given
L3 = Mv*deltav;



%2nd set up grid
%grip is in 3 dimensions, v-dimension, s-dimension, tau-dimension

deltas = deltav * tau;
smax = si(end); %the largest s element
Ms = smax/deltas; %elements of s-grid
%Mv = elements of v-grid is already given
Mtau = ...; 
  
%todo setup grid, dynamic of fixed length?
c = ...;

%3rd add costs for grid
v[k+1] = v[k] + 1/deltav *(c[1,k]*tau + c[2,k] *tau*tau/2);
s[k+1] = s[k] + v[k] + 1/deltas *(c[1,k]*tau*tau/2 + c[2,k] *tau*tau*tau/6);


%todo solve to c
c2[k] = (s[k+1]-s[k]-v[k])*deltas - (v[k+1]-v[k])*deltav*tau/2)/(tau*tau*tau*(1/6-1/4));
c1[k] = (v[k+1]-v[k])*deltav/tau-c2[k]*tau/2;
%use limits
%4th find optimum path



end

