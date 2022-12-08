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
    mul_2 = transpose(r{2}(id_TakeDecision(:,1)+ (id_TakeDecision(:,2)-1)*n))

    % get current transition
    factor_multiply=mul_1.*(sum_id(:, [1, end])==0)+mul_2.*(sum_id(:, [1, end])==1)+(sum_id(:, [1, end])==2);
  
    P_Prod_Check=ones(size(P(P(:, i)==1, i), 1), 1); % Make sure its valid
    
    % if N greater then N
    if N>2
        id_TakeDecision=[(MatStateID(i, 1:end-2)+1)', (MatStateID(i, 3:end)+1)']; % end state
        mul1 = transpose(r{1}(id_TakeDecision(:,1) + (id_TakeDecision(:,2)-1)*n));
        mul2 = transpose(r{2}(id_TakeDecision(:,1)+ (id_TakeDecision(:,2)-1)*n))
        
        P_elements=mul1.*(sum_id(:, 2:end-1)==0)+mul2.*(sum_id(:, 2:end-1)==1)+(sum_id(:, 2:end-1)==2);
        for k=1:size(P_elements, 2)
            P_Prod_Check=P_Prod_Check.*P_elements(:, k); % transition
        end
    end
    P(P(:, i)==1, i)=P_Prod_Check.*factor_multiply(:, 1).*factor_multiply(:, 2);
end

cost = ones(hbar+1, 2^N)*delta; % start with a full of grid which wants to terminate constantly
cost_prev = ones(hbar, 2^N)*delta; % start with a full of grid which wants to terminate constantly

policy = zeros(hbar, TotalStateAmount);
beta_u =75;
beta_d = 0;

beta = beta_u;
corr = false;
while ( ~corr)
    beta = (beta_u+beta_d)/2;
for h_ =(hbar):-1:1
    for i =1:(TotalStateAmount)
      cur = MatStateID(i,:);
     sun_ones_curr_state =sum(MatStateID(i,:),2)*d;
     P_curr =P(:,i);
     cut_cost = cost(h_+1,:) ;
   %  cur_cost_prev = cost_prev(i,:);
   cost_j = dot(P_curr,cut_cost);
     s_cost = dot(cut_cost,(P_curr)) + sun_ones_curr_state - beta;
    cost(h_,i) =  min([s_cost,delta]);
    end
end
    if(cost(1,1)>0)
    beta_d = beta;
    else
      beta_u = beta;
    end
    if(abs(cost(1,1))<=0.0001)
    corr = true;
    end
end


ssiu =1;

for i=1:hbar+1
    for j = 1 :TotalStateAmount
    if(cost(i,j)>= delta)
        cost(i,j) = delta; 
        policy(i,j) = 1;
    end
    end
end

matrix_output = policy;

for i =1:hbar+1
matrix_output(i,:) = matrix_output(i,:)*i;
end
for h_ = (hbar):-1:1
for next =1:TotalStateAmount
     P_curr =P(:,next);
     cur_mat_value = matrix_output(h_, next);
     if(matrix_output(h_, next)==0)
      next_policy = matrix_output(h_+1, :);
    matrix_output(h_, next) = dot(P_curr, next_policy);
     else
         matrix_output(h_, next) = h_;
     end
end
end
Js = beta;
avtau = matrix_output(1,1);
mu_ = policy(2:end-1,:)'+1;