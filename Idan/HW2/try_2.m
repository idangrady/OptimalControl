clear

clear

delta = 7;
d     = 1;
r     = [0.1 0.2 0.6];
hbar  = 10;
N     = 2;
ep = 0.1; ell = 1; gamma = 1;

TotalStateAmount = 2^N;
r_1_1=r(3); r_0_0 =r(1);r_0_1 = r(2);

r={[1-r_0_0, 1-r_0_1; 1-r_0_1, 1-r_1_1], [r_0_0, r_0_1; r_0_1, r_1_1]};
n=size(r{1}, 1); 
id=0:2^N-1; 

idBinary=dec2bin(id); 
MatStateID=[];

for i=1:N
    MatStateID=[MatStateID, str2num(idBinary(:, i))]; 
end
P=ones(size(idBinary, 1), size(idBinary, 1)); 
for i=1:size(idBinary, 1) 
    current = MatStateID(i, :);
    all = MatStateID;
    current_all = all - current;
    
    % make sure 1-> 0 transition is impossible
    sum_cur_all = current_all==-1;
    cur_sum =sum(sum_cur_all, 2);
    index = find(cur_sum==1);
    if(~isempty(index))
        P(index, i)=0;
    end

    
    wer = P(:, i)==1;
    sum_id=MatStateID(P(:, i)==1, :)+MatStateID(i, :);

    id_TakeDecision=[[1;1], (MatStateID(i, [2, end-1])+1)']; 
    
    mul_1 = transpose(r{1}(id_TakeDecision(:,1) + (id_TakeDecision(:,2)-1)*n)); % from [0 to 1]
    mul_2 = transpose(r{2}(id_TakeDecision(:,1)+ (id_TakeDecision(:,2)-1)*n));

    % get current transition
    factor_multiply=mul_1.*(sum_id(:, [1, end])==0)+mul_2.*(sum_id(:, [1, end])==1)+(sum_id(:, [1, end])==2);
  
    P_Prod_Check=ones(size(P(P(:, i)==1, i), 1), 1); % Make sure its valid
    
    % if N greater then N
    if N>2
        id_TakeDecision=[(MatStateID(i, 1:end-2)+1)', (MatStateID(i, 3:end)+1)']; % end state
        mul1 = transpose(r{1}(id_TakeDecision(:,1) + (id_TakeDecision(:,2)-1)*n));
        mul2 = transpose(r{2}(id_TakeDecision(:,1)+ (id_TakeDecision(:,2)-1)*n));
        
        P_elements=mul1.*(sum_id(:, 2:end-1)==0)+mul2.*(sum_id(:, 2:end-1)==1)+(sum_id(:, 2:end-1)==2);
        for k=1:size(P_elements, 2)
            P_Prod_Check=P_Prod_Check.*P_elements(:, k); % transition
        end
    end
    P(P(:, i)==1, i)=P_Prod_Check.*factor_multiply(:, 1).*factor_multiply(:, 2);
end


% Define the cost-to-go matrix C, which describes the costs associated
% with transitioning between states
C = ones(hbar, TotalStateAmount);

% Define the sensor measurement error probability ep

% Initialize the prior probabilities for each state
prior = ones(1, TotalStateAmount)/ TotalStateAmount;

% Define the observation likelihood matrix L, which describes the
% probability of observing a particular measurement given each state
L = [ones(1, TotalStateAmount)*ep; ones(1, TotalStateAmount)*(1-ep)]; % the measurement of the sensor
policy = ones(hbar, TotalStateAmount);
for i=1: hbar
% Define the initial observation to be made
initial_observation = 1;
id_1s = [0 1 1 2];

for state=1:TotalStateAmount
% Calculate the likelihood of the initial observation given each state


next_state_probability = P(:,i);
likelihood = L(:, initial_observation+1);
next_transition_cost = 

% Use Bayes' theorem to update the posterior probabilities for each
% state, using the initial likelihood and the prior probabilities
posterior = prior.* transpose(next_state_probability);
posterior = posterior / sum(posterior); % divide by alpha

% Use the transition matrix P and the cost-
  % the expected cost-to-go for each state, assuming that the current
  % state is the one with the highest posterior probability
  [~, max_idx] = max(posterior);
  expected_cost = C(i,:) * posterior;

  % Update the cost-to-go matrix C based on the observed transitions
  % and the resulting expected costs-to-go
  C = update_costs(C, max_idx, expected_cost);

  % Update the posterior probabilities based on the updated cost-to-go
  % matrix C
  posterior = expected_cost .* posterior;
  posterior = posterior / sum(posterior);
end

end

% The final posterior probabilities give the estimated expected cost-to-go
% for each state given the sequence of observations
disp(posterior);

% Define a function for updating the cost-to-go matrix based on the
% observed transitions and the resulting expected costs-to-go
function C = update_costs(C, max_idx, expected_cost)
  % Update the costs associated with transitioning to the state with
  % the highest posterior probability, based on the expected cost-to-go
  % for that state
  C(:, max_idx) = expected_cost;
end
