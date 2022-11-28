function [c] = timeoptimalgridCheck(omegai,si,tau,deltav,L2,Na,Mv)

L1 = deltav*Na/tau;
L3 = Mv*deltav;
deltas = deltav * tau;
sl = si(end); 
Ms = sl/deltas; 
h = length(si);   

Mv = 2;
Ms = 2;
block =(1+Mv)*(Ms+1);
grid = zeros((Mv+1)*(Ms+1),(1+Ms)*(1+Mv));
C_1 = zeros((Mv+1)*(Ms+1),(1+Ms)*(1+Mv));
C_2 = zeros((Mv+1)*(Ms+1),(1+Ms)*(1+Mv));

grid_check = zeros(Mv+1,Ms+1);
grid_next_check = zeros(Mv+1, Ms+1);
for curr_state =1:block % assumed Mv comes first. 
for next_state = 1:block

     [v_current, s_current] = idx_to_state_velo(curr_state, Mv, Ms);
     [v_next,s_next] = idx_to_state_velo(next_state, Mv, Ms);
      s_next = s_next+1; s_current = s_current+1;
        grid_check(v_current, s_current) = curr_state;
        grid_next_check(v_next,s_next) = next_state;
    

end
end


end