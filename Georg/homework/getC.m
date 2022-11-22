function [c] = getC(v_one,s_one,tau,deltav)

c = zeros(2,length(v_one));
deltas = deltav * tau;
for i = 1:length(v_one)

  if(i<length(v_one))
    v_current = v_one(i);
    v_next = v_one(i+1);
    s_current = s_one(i);
    s_next = s_one(i+1);

    c(2,i) = ((s_next-s_current-v_current)*deltas - (v_next-v_current)*deltav*tau/2)/(tau^3*(1/6-1/4));
    c(1,i) =(v_next-v_current)*deltav/tau-c(2,i)*tau/2;
    
  else
    c(1,i) = c(1,i-1);
    c(2,i) = c(2,i-1);
  end
      
end

