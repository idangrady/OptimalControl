m = 5; % block-> in our case is Mv*Ms

a = [0 0 inf inf 1; 3 0 5 3 1; inf 5 0 1 3 ; inf inf 1 0 1 ; 1 1 3 1 0 ]; % grid
h= 5; % maximm depth

end_state = 3;  % to fill from function
start_state = 1; % to fill in from function
% start from end
cost = a(:, end_state);
policy = inf(1,5);

for i= 1:m % node m = Mv*Ms
    if(a(i,end_state) == inf)
        continue;
    else
        policy(i) = end_state;
    end
end

for cur_ =1:m
    cur_state_space = a(:,cur_);
    p = length(cur_state_space);
    for previous = 1:length(cur_state_space)
    if(cost(cur_) + cur_state_space(previous)< cost(previous))
        policy(previous) = cur_;
        cost(previous) =cost(cur_) + cur_state_space(previous);
    end
    end
end


% traverse the path
output = [start_state];
next_step=start_state;
for i=1:h
    
next_step = policy(next_step);
output(end+1) = next_step; % output is the currect sequence of steps
if(next_step ==end_state)
break
end
end

% traverse the path
% need to convert back the idx state into velocity and state