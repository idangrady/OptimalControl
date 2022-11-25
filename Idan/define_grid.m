function [C,M,number_elements] = define_grid(omegai,si,tau,deltav,L2,Na,Mv,enable_constraints)

enable_constraints=true;
L1 = deltav*Na/tau;
L3 = Mv*deltav;
deltas = deltav * tau;
sl = si(end); 
Ms = sl/deltas +1; 
h = length(si);   
number_elements = Mv*Ms;
grid = zeros(Mv*Ms,Ms*Mv);


%problem with Ms
v = linspace(0, Mv, Mv+1) ; % I changed the multiplication as I dont believe it is necessary
s = linspace(0, Ms, Ms+1);  % I changed the multiplication as I dont believe it is necessary
i =0;
o =0;
for v_currentIdx = 1:Mv 
  for s_current_Idx = 1:Ms
      i=i+1;
      o =1;
    for v_nextIdx = 1:Mv 
      for s_nextIdx = 1:Ms
            
        s_current = s(s_current_Idx);
        s_next    = s(s_nextIdx);
        v_current = v(v_currentIdx);
        v_next    = v(v_nextIdx);

        %todo - think about initial state

            c2 = ((s_next-s_current-v_current)*deltas - (v_next-v_current)*deltav*tau/2)/(tau^3*(1/6-1/4));
            c1 = (v_next-v_current)*deltav/tau-c2*tau/2;

            if((c2<0) && (c1>0))
              vintmax = v_current-c1^2/2/c2;
            else
              vintmax = v_current;
            end
            
            if(enable_constraints)
              if(abs(c1)>L1)
                cost = inf;
              elseif(abs(c1+c2*tau)>L1)
                cost = inf;
              elseif(omegai(k)*vintmax^2 > L2)
                cost = inf;
              elseif(vintmax>L3)
                cost = inf;
              else
    %             cost = c1 + c2;   %maybe wrong way to define the cost
                cost =1; % (si(k) - s_current)^2; %assumption about quadratic error
              end
            else 
               cost = 1;% (si(k) - s_current)^2; %assumption about quadratic error
            end
    
            
            %beginning contraint
            if(k==1) %s(0)=0 && %v(0)=0
              if(not((s_current==0) && (v_current==0)))
                cost = inf;
              end
            end
            
            %end contraint
            if(k==h) %s(h)=sl %v(0)=0
              if(not(v_current==0))
                cost = inf;
              end
            end
            
            %cannot go backwards
            if(s_next<s_current)
              cost = inf;
            end

        %assign cost to grid
        grid(i,o) = cost;

      %  grid(s_nextIdx+(v_nextIdx-1)*Mv+(s_current_Idx-1)*(Ms*Mv)+ (v_currentIdx-1) * (Ms*Mv*Ms),k) = cost;
       o=o+1;

      end
    end
  end
end

start_idx =2;
goal_idx = Mv*Ms = 2; % this is not hard coded yet will be taken from the mv and Ms state;

end_space = grid(:,goal_idx);
cost_func =grid(:,goal_idx);

policy = zeros(1,Ms*Mv);
AccumelatedCost=(1,Mv*Ms);

max_k = 10;
for k 1:h % 
  for currState = 1:Ms*Mv
       if(start_idx ==goal_idx)
       % break
       return_val = 0; % change that to fit the constraint
       break
       end

      % check whether we exceeded the max itteration and there is no path
      % if K> maxk return false

      % state space to get to Last node
      exploreSpace = grid(:,currState);
      for i 1:length(exploreSpace)
          if(cost_func(i)==inf) % more efficient as it can not be updated
            continue;
          end
      % each i corresponds to a specific index of state, velocity space. 
          
          if(grid(i,currState)) ==inf
              continue;
      end
      % go from the end to the beginning. 

  
  end
end

end

end
