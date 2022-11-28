function [c] = timeoptimalpathspeed2(omegai,si,tau,deltav,L2,Na,Mv)
%parameters:
%omegai     vector,   given omegavector
%si...      vector,   given distance vector
%tau...     constant, stepsize for the time
%deltav...  constant, stepsize for speed
%L2...      constant, limit for velocity change
%Na...      constant
%Mv...      constant

%% first: set up grid
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
  h = length(si);   %expected that this is possible

  %problem with Ms
  v = 0:deltav:Mv*deltav;
  s = 0:deltas:Ms*deltas;
  s_normed = 0:1:Ms;
  v_normed = 0:1:Mv;
  C = ones(h,(Mv+1)*(Ms+1),(Mv+1)*(Ms+1))*inf;
  M = zeros(h,(Mv+1)*(Ms+1),(Mv+1)*(Ms+1));

  %fill grid
  for k = 1:h
    for v_currentIdx = 1:Mv+1 
      for s_current_Idx = 1:Ms+1
        for v_nextIdx = 1:Mv+1 
          for s_nextIdx = 1:Ms+1
            s_normed_current = s_normed(s_current_Idx);
            s_normed_next    = s_normed(s_nextIdx);
            v_normed_current = v_normed(v_currentIdx);
            v_normed_next    = v_normed(v_nextIdx);
            c2 = ((s_normed_next-s_normed_current-v_normed_current)*deltas - (v_normed_next-v_normed_current)*deltav*tau/2)/(tau^3*(1/6-1/4));
            c1 = (v_normed_next-v_normed_current)*deltav/tau-c2*tau/2;

            if((c2<0) && (c1>0))
              vintmax = v_normed_current-c1^2/2/c2;
            else
              vintmax = v_normed_current;
            end
          
            if(abs(c1)>L1)
              cost = inf;
            elseif(abs(c1+c2*tau)>L1)
              cost = inf;
            elseif(omegai(k)*vintmax^2 > L2)
              cost = inf;
            elseif(vintmax>L3)
              cost = inf;
            else
              cost = (si(k) - s_normed_current*deltas)^2; %assumption about quadratic error
            end
   
            %beginning contraint
            if(k==1) %s(0)=0 && %v(0)=0
              if(not((s_normed_current==0) && (v_normed_current==0)))
                cost = inf;
              end
            end
            
            %end contraint
            if(k==h) %s(h)=sl %v(0)=0
              if(not(v_normed_current==0))
                cost = inf;
              end
            end
          
            %cannot go backwards
            if(s_normed_next<s_normed_current)
              cost = inf;
            end
           
            %assign cost to grid
            C(k,s_current_Idx+Ms*(v_currentIdx-1),s_nextIdx+Ms*(v_nextIdx-1)) = cost;
            M(k,s_current_Idx+Ms*(v_currentIdx-1),s_nextIdx+Ms*(v_nextIdx-1)) = s_nextIdx+Ms*(v_nextIdx-1);          
          end
        end
      end
    end
  end
end


