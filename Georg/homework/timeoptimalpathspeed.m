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
  sl = si(end); %the largest s element
  Ms = sl/deltas; %elements of s-grid
  %Mv = elements of v-grid is already given

  %problem with Ms
  v = 0:deltav:Mv*deltav;
  s = 0:deltas:Ms*deltas;
  s_normed = 0:1:Ms;
  v_normed = 0:1:Mv;
  C1 = ones(Ms+1,Mv+1)*inf;
  C2 = ones(Ms+1,Mv+1)*inf;
  transition_possible = ones(Ms+1,Mv+1);
  v_tmp =ones(Ms+1,Mv+1);
  s_chosen = [0];
  v_chosen = [0];
  c = zeros(2,1);
  h = 0;
  
  s_normed_current = 0;
  v_normed_current   = 0;
  skip_loop = 0;
  i = 1;
 
  diff = abs(sl - s_normed_current*deltas)
  while( diff >deltas)   %repeat until goal is reached
    diff = abs(sl - s_normed_current*deltas)
    h = h+1;
    if((s_normed_current*deltas >= si(i)) && (i<length(omegai)))
      i = i+1;
    end
    for s_kp1 =s_normed_current+1:Ms+1
      for v_kp1 = 1:Mv+1
         C2(s_kp1,v_kp1) = ((s_kp1-s_normed_current-v_normed_current)*deltas - (v_kp1-v_normed_current)*deltav*tau/2)/(tau^3*(1/6-1/4));
         C1(s_kp1,v_kp1) = (v_kp1-v_normed_current)*deltav/tau-C2(s_kp1,v_kp1)*tau/2;
%            C1_tmp = ((s_kp1-s_normed_current-v_normed_current)*deltas - (v_kp1-v_normed_current)*tau/3*deltav)/tau^2*6;
         %check if transition is possible
         if((C2(s_kp1,v_kp1)<0) && (C1(s_kp1,v_kp1)>0))
            vintmax = v(v_normed_current+1)-C1(s_kp1,v_kp1)^2/2/C2(s_kp1,v_kp1);
          else
            vintmax = v(v_normed_current+1);
         end
         v_tmp(s_kp1,v_kp1) = vintmax;
          if(abs(C1(s_kp1,v_kp1))>L1)
            transition_possible(s_kp1,v_kp1) = 0;
          end
          if(abs(C1(s_kp1,v_kp1)+C2(s_kp1,v_kp1)*tau)>L1)
            transition_possible(s_kp1,v_kp1) = 0;
          end
          if(omegai(i)*vintmax^2 > L2)
            transition_possible(s_kp1,v_kp1) = 0;
          end
          if(vintmax>L3)
            transition_possible(s_kp1,v_kp1) = 0;
          end
      %every possible transition will be remaining with an 1
      end  
    end
    %choose optimal state
    for s_kp1 =Ms:-1:1  %try to make as big steps as possible
      for v_kp1 = 1:Mv+1      %try to have as little speed as possible on the check points
        if(transition_possible(s_kp1,v_kp1) ==1)
          s_chosen(h) = s_kp1;
          v_chosen(h) = v_kp1;
          s_normed_current = s_kp1;
          v_normed_current = v_kp1;
          c(1,h) = C1(s_kp1,v_kp1);
          c(2,h) = C2(s_kp1,v_kp1);
          skip_loop =1;
          break;
        end
      end
      if(skip_loop==1)
        break;
      end
    end %vkp1 for
    skip_loop = 0;
    transition_possible = ones(Ms+1,Mv+1);
    v_tmp =ones(Ms+1,Mv+1);

  end %while


  c(1,h+1) = c(1,h);
  c(2,h+1) = c(2,h);

  hold on
  plot(v_chosen,"g")
  plot(s_chosen,"b-")
%   plot(,"rx")
  grid on;
  hold off
  legend("speed","optimal distance","wished distance")
  
end

