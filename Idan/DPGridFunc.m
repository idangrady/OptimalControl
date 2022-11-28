function [output, policy, cost] = DPGridFunc(grid, end_state, start_state, max_depth) % V_imax ,L2
%grid = [0 0 inf inf 1; 3 0 5 3 1 ; inf 5 0 1 3 ; inf inf 1 0 1; 1 1 3 1 0];

m = length(grid); % block-> in our case is Mv*Ms

cost = grid(:, end_state);
policy = inf(1,m);

for i= 1:m % node m = Mv*Ms
    if(grid(i,end_state) == inf)
        continue;
    else
        policy(i) = end_state;
    end
end

for cur_ =m:-1:1
    cur_state_space = grid(:,cur_);
    p = length(cur_state_space);
    for previous = 1:length(cur_state_space)
    if(cost(cur_) + cur_state_space(previous)< cost(previous))
        policy(previous) = cur_;
        cost(previous) =cost(cur_) + cur_state_space(previous);
    end
    end
end

% traverse the path
output = [start_state]; %start_state
for i=1:max_depth
next_step = policy(output(i));
output(end+1) = next_step; % output is the currect sequence of steps

if(next_step ==end_state)
break

end
end

end