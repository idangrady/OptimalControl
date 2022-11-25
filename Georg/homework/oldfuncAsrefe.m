% here come the function//

enable_constraints=true;
L1 = deltav*Na/tau;
L3 = Mv*deltav;
deltas = deltav * tau;
sl = si(end); 
Ms = sl/deltas +1; 
h = length(si);   
number_elements = Mv*Ms;
grid = zeros(Mv*Ms,Ms*Mv);

Mgrid = zeros(Mv*Ms*Ms*Mv, h);

%problem with Ms
v = linspace(0, Mv, Mv+1) ; % I changed the multiplication as I dont believe it is necessary
s = linspace(0, Ms, Ms+1);  % I changed the multiplication as I dont believe it is necessary

for k = 1:h
i =1;
for v_currentIdx = 1:Mv 
  for s_current_Idx = 1:Ms
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
       Cgrid(i,k) = cost;
       Mgrid(i,k) = s_nextIdx+Ms*(v_nextIdx-1);
        i=i+1;
      %  grid(s_nextIdx+(v_nextIdx-1)*Mv+(s_current_Idx-1)*(Ms*Mv)+ (v_currentIdx-1) * (Ms*Mv*Ms),k) = cost;
      end
    end
  end
end
end


block = Ms*Mv;

M_best_path_idx = zeros(Ms*Mv,h);
policyM = zeros((Mv*Ms),h); 

policyV = zeros((Mv*Ms),h); 
V_best_path_idx = zeros(Ms*Mv,h);

action= zeros((Mv*Ms),h);
% Loop over all the next state and compare and find the min value between
% them. 
for k = 0:h-1 % go over all the maximu, st

    cur_k =cast(h-k,"uint8"); % lol not sure how to do it like in python 

    for currstate= 1:Mv*Ms
    % Every block of Mv*Mv we only change the input variable ((0,0),(1,0)) ->
    % ((0,1),(1,0_) only the origin state changes. Hence we need to find
    % the smallest between those
    
    check_next_state = cast(currstate+ (currstate-1)*block,"uint8");
    MS_val = Mgrid(check_next_state, cur_k);
    if(policyM(check_next_state, cur_k)==0)
        policyM(check_next_state, cur_k) = MS_val;
        M_best_path_idx(currstate,cur_k) = (currstate-1)*block; % gives the ouput idx direct

    elseif(MS_val < policyM(check_next_state, cur_k))
        policyM(check_next_state, cur_k) = MS_val;
        M_best_path_idx(currstate,k) = (currstate-1)*block; % gives the ouput idx direct
    end
end

end

end
