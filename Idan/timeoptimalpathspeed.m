function [c] = timeoptimalpathspeed(omegai,si,tau,deltav,L2,Na,Mv)

L1 = deltav*Na/tau;
L3 = Mv*deltav;
deltas = deltav * tau;
sl = si(end); 
Ms = sl/deltas; 
h = length(si);   

grid = zeros((Mv+1)*(Ms+1),(1+Ms)*(1+Mv));
C_1 = zeros((Mv+1)*(Ms+1),(1+Ms)*(1+Mv));
C_2 = zeros((Mv+1)*(Ms+1),(1+Ms)*(1+Mv));

V_imax = zeros((Mv+1)*(Ms+1),(1+Ms)*(1+Mv));
block = (Mv+1)*(Ms+1);

last_state = block -(Ms);
begin_state = 1;

i =1;
for curr_state =1:block % assumed Mv comes first. 
for next_state = 1:block

     [v_current, s_current] = idx_to_state_velo(curr_state, Mv);
     [v_next,s_next] = idx_to_state_velo(next_state, Mv);
     
     s_current = deltas*(s_current);s_next =deltas*(s_next); % try this
     v_next = deltav*(v_next);v_current =deltav*(v_current); % try this
     
     c2 = ((s_next-s_current-v_current)*deltas - (v_next-v_current)*deltav*tau/2)/(tau^3*(1/6-1/4));
     c1 = (v_next-v_current)*deltav/tau-c2*tau/2;

      if((s_current*deltas >= si(i)) && (i<length(omegai)))
      i = i+1;
      end

     %for future reference
     C_1(curr_state, next_state) = c1;
     C_2(curr_state, next_state) = c2;

            if((c2<0) && (c1>0))
              vintmax = v_current-c1^2/2/c2;
            else
              vintmax = v_current;
            end

              if(abs(c1)>L1)
                cost = inf;
              elseif(abs(c1+c2*tau)>L1)
                cost = inf;
              elseif(omegai(i)*vintmax^2 > L2)   
                cost = inf;
              elseif(vintmax>L3)
                cost = inf;
              elseif (curr_state==next_state)
            cost = inf;
              else
                cost =1; 
              end
          %  if(cost<inf)
          %  fprintf("Now");
          %  end
        grid(curr_state,next_state) = cost;

end
end

enState = length(grid)-Ms;
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