function [c] = reference(omegai,si,tau,deltav,L2,Na,Mv)

L1 = deltav*Na/tau;
L3 = Mv*deltav;
deltas = deltav * tau;
sl = si(end); 
Ms = sl/deltas +1; 
h = length(si);   

grid = zeros((Mv+1)*(Ms+1),(1+Ms)*(1+Mv));
C_1 = zeros((Mv+1)*(Ms+1),(1+Ms)*(1+Mv));
C_2 = zeros((Mv+1)*(Ms+1),(1+Ms)*(1+Mv));

%
c = 0;
V_imax = zeros((Mv+1)*(Ms+1),(1+Ms)*(1+Mv));
block = (Mv+1)*(Ms+1),(1+Ms)*(1+Mv);
last_state = ;%si(end);
begin_state = 3;%si(1);

for curr_state =1:block % assumed Mv comes first. 
for next_state = 1:block

     [v_current, s_current] = idx_to_state_velo(curr_state, Mv, Ms);
     [v_next,s_next] = idx_to_state_velo(next_state, Mv, Ms);

     c2 = ((s_next-s_current-v_current)*deltas - (v_next-v_current)*deltav*tau/2)/(tau^3*(1/6-1/4));
     c1 = (v_next-v_current)*deltav/tau-c2*tau/2;
     
     %for future reference
     C_1(curr_state, next_state) = c1;
     C_2(curr_state, next_state) = c2;

            if((c2<0) && (c1>0))
              vintmax = v_current-c1^2/2/c2;
            else
              vintmax = v_current;
            end

            %V_imax(curr_state, next_state) =vintmax;
                % constraints
              if(abs(c1)>L1)
                cost = inf;
              elseif(abs(c1+c2*tau)>L1)
                cost = inf;
%              elseif(omegai(k)*vintmax^2 > L2)   How to bring k back?
%                cost = inf;
              elseif(vintmax>L3)
                cost = inf;
              else
                cost =1; 
              end
    
            if (curr_state==next_state)
            cost = 0;
            end
            
            % an idea is to assign that 
        %assign cost to grid
        grid(curr_state,next_state) = cost;

end
end

enState = length(grid);
begState = 1;
c =0;
[path, policy, cost] = DPGridFunc(grid,enState, begState, h ); % V_imax, L2
outc_1 = [];
outc_2 = [];
for i  =1:length(path)
outc_1(end +1) =C_1(path(i), path(i+1));
outc_2(end +1) =C_2(path(i), path(i+1));
end

while(length(outc_1)<h)
    outc_1(end+1) = outc_1(end);
    outc_2(end+1) = outc_2(end);
end

end