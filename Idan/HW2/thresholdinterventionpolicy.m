function [Js,avtau]= thresholdinterventionpolicy(delta,d,r,hbar,N,ep,ell,gamma)

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
    index = find(cur_sum>=1);
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

permut = 0:2^hbar-1;
Permutation=dec2bin(permut); 

% compute possible permutation
permut_Mat=[];

for i=1:(hbar)
    permut_Mat=[permut_Mat, str2num(Permutation(:, i))]; 
end


% Define the cost-to-go matrix C, which describes the costs associated
% with transitioning between states
C = zeros( 2^hbar, hbar+1);
policy = zeros( 2^hbar, hbar);
R = zeros(2,2^N);

%assume we could simply take the ell measurement
for i=1:TotalStateAmount
    cur = MatStateID(i, ell);
    if(cur ==1)
        R(2,i)=(1-ep);
        R(1,i)=(ep);
    else
      R(1,i)=(1-ep);
      R(2,i)=(ep);
    end
end


p_klein =ones(2^hbar, hbar);
d_mat = sum(MatStateID, 2)'.*d;

% bayes filter
for s =1: 2^hbar
prior =zeros(TotalStateAmount,1);prior(1)=1; % start 
intervene = false;
for h_=1: hbar
    if(h_>1)
    p_k =P*prior; % estimation step
    else 
        p_k = prior;
    end
    given_val = permut_Mat(s, h_)+1;

    % get the sensor measurement
    r_ = (R(given_val, :).*eye(TotalStateAmount));

    %update step
    sum_probability = sum(r_*p_k);
    if(h_>1)
    p_klein(s, h_) =p_klein(s, h_-1)* sum_probability;
    else
     p_klein(s, h_) = sum_probability;
    end
    prior = r_*p_k / sum_probability; %
    if(~intervene)
    C(s, h_) = d_mat* prior;

    else
       C(s, h_) = 0;
    end
    if(C(s, h_)>gamma)
     C(s, h_) =delta;
     policy(s,  h_-1)=1;

     intervene = true;
    end
end

if(~intervene)
    C(s, h_+1)=delta;
    policy(s,  h_)=1;
end
end 


matrix_output = policy;
p_klein = p_klein./sum(p_klein);

arr = 1:hbar;
matrix_output = matrix_output*arr';
last_colPklein = p_klein(:,end);
avtau = sum(last_colPklein.*matrix_output);
Js = dot(sum(C, 2)',last_colPklein);

end