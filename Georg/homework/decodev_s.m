function [v_one,s_one] = decodev_s(si,tau,deltav,Mv,u_one)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

  deltas = deltav * tau;
  sl = si(end); %the largest s element
  Ms = sl/deltas +1; %elements of s-grid
  v = linspace(0, (Mv-1)*deltav, Mv) ;
  s = linspace(0, (Ms-1)*deltas, Ms);
  len = length(u_one);
  s_one = zeros(len,1);
  v_one = zeros(len,1);
  
  s_one(1) = 0;
  v_one(1) = 0;
  
  for i = 2:len
    idx_v = floor((u_one(i-1)-1)/Ms)+1;
    idx_s = u_one(i-1)-Ms*(idx_v-1);
    v_one(i) = v(idx_v);
    s_one(i) = s(idx_s);
  end
  
  
end

