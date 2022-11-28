function [c] = timeoptimalpathspeednew(omegai,si,tau,deltav,L2,Na,Mv)


L1 = deltav*Na/tau;
L3 = Mv*deltav;
deltas = deltav * tau;
sl = si(end); 
Ms = sl/deltas +1; 
h = length(si);   

number_elements = Mv*Ms;

grid = zeros((Mv+1)*(Ms+1),(1+Ms)*(1+Mv));
C_1 = zeros((Mv+1)*(Ms+1),(1+Ms)*(1+Mv));
C_2 = zeros((Mv+1)*(Ms+1),(1+Ms)*(1+Mv));

last_state = si[end];
begin_state = si[1];


m = length(grid); % block-> in our case is Mv*Ms

cost = grid(:, end_state);
policy = inf(1,5);

for i= 1:m % node m = Mv*Ms
    if(grid(i,end_state) == inf)
        continue;
    else
        policy(i) = end_state;
    end
end

for curr_state =(Mv+1)*(Ms+1):-1:1 % assumed Mv comes first.  h:-1:1
        
    cur_state_space = grid(:,curr_state);
  %  for previous = 1:length(cur_state_space)
    for next_state = 1:length(cur_state_space)

     [v_current, s_current] = idx_to_state_velo(curr_state, Mv, Ms);
     [v_next,s_next] = idx_to_state_velo(next_state, Mv, Ms);

     c2 = ((s_next-s_current-v_current)*deltas - (v_next-v_current)*deltav*tau/2)/(tau^3*(1/6-1/4));
     c1 = (v_next-v_current)*deltav/tau-c2*tau/2;


    

    if(cost(cur_) + cur_state_space(previous)< cost(previous))
        policy(previous) = cur_;
        cost(previous) =cost(cur_) + cur_state_space(previous);
    end
    end
end


% traverse the path
output = [start_state];
next_step=start_state;
for i=1:max_depth
    
next_step = policy(next_step);
output(end+1) = next_step; % output is the currect sequence of steps
if(next_step ==end_state)
break
end
end



grid = zeros((Mv+1)*(Ms+1),(1+Ms)*(1+Mv));
C_1 = zeros((Mv+1)*(Ms+1),(1+Ms)*(1+Mv));
C_2 = zeros((Mv+1)*(Ms+1),(1+Ms)*(1+Mv));


end;