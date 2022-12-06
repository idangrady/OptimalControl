
clc()
clear

delta = 7;
d     = 1;
r     = [0.1 0.2 0.6];
hbar  = 10;
N     = 2;
ep = 0.1; ell = 1; gamma = 1;


r_11=r(3); r_00 =r(1);r_01 = r(2);
TotalStateAmount = 2^N;
r={[1-r_00, 1-r_01; 1-r_01, 1-r_11], [r_00, r_01; r_01, r_11]};
id=0:2^N-1; %index numbers
n=size(r{1}, 1); %size of r
idBinary=dec2bin(id); % Convert index number in decimal integer to its binary representation
MatStateID=[];



for i=1:N
    MatStateID=[MatStateID, str2num(idBinary(:, i))]; % Binary representation of the index number in matrix form
end
P=ones(size(idBinary, 1), size(idBinary, 1)); % P matrix; Size of the matrix is determined by 2^N
for i=1:size(idBinary, 1) % size of rows
    current = MatStateID(i, :);
    all = MatStateID;
    current_all = current-all;

    
    P(sum((MatStateID-MatStateID(i, :))==-1, 2)>0, i)=0; %Substitute P{x_{t+1, l}=0 | x_{t, l}=1}=0

    wer = P(:, i)==1;
    add_id=MatStateID(P(:, i)==1, :)+MatStateID(i, :); %Find transitions between stages(add_id=1 when 0->1, add_id=0 when 0->0 and add_id=2 when 1->1)
    id_TakeDecision=[[1;1], (MatStateID(i, [2, end-1])+1)']; %Set of index numbers to call r for probability computation when P{x_{t+1, 1}=1 | x_{t, 1}=0, x_{t, 2}=a} and P{x_{t+1, N}=1 | x_{t, N-1}=a, x_{t, N}=0}, etc.
    PE_Multiply_01=transpose(r{1}(id_TakeDecision(:,1) + (id_TakeDecision(:,2)-1)*n)).*(add_id(:, [1, end])==0)+transpose(r{2}(id_TakeDecision(:,1)...
        + (id_TakeDecision(:,2)-1)*n)).*(add_id(:, [1, end])==1)+(add_id(:, [1, end])==2);%Probability computation when P{x_{t+1, 1}=1 | x_{t, 1}=0, x_{t, 2}=a} and P{x_{t+1, N}=1 | x_{t, N-1}=a, x_{t, N}=0}, etc.
    
    P_Prod_Check=ones(size(P(P(:, i)==1, i), 1), 1);
    if N>2
        id_TakeDecision=[(MatStateID(i, 1:end-2)+1)', (MatStateID(i, 3:end)+1)']; %Set of index numbers to call r for probability computation when P{x_{t+1, l}=1 | x_{t, l-1}=a, x_{t, l}=0, x_{t, l+1}=b}
        P_elements=transpose(r{1}(id_TakeDecision(:,1) + (id_TakeDecision(:,2)-1)*n)).*(add_id(:, 2:end-1)==0)+transpose(r{2}(id_TakeDecision(:,1)...
            + (id_TakeDecision(:,2)-1)*n)).*(add_id(:, 2:end-1)==1)+(add_id(:, 2:end-1)==2); %Probability computation when P{x_{t+1, l}=1 | x_{t, l-1}=a, x_{t, l}=0, x_{t, l+1}=b}
        for k=1:size(P_elements, 2)
            P_Prod_Check=P_Prod_Check.*P_elements(:, k);
        end
    end
    P(P(:, i)==1, i)=P_Prod_Check.*PE_Multiply_01(:, 1).*PE_Multiply_01(:, 2);%Combine all the probability
end


cost = ones(1,2^N)*delta; % start with a full of grid which wants to terminate constantly
cost_prev = ones(1,2^N)*delta; % start with a full of grid which wants to terminate constantly

policy = zeros(hbar, TotalStateAmount);

for h_ =1:hbar
    for i =1:length(cost)
     sun_ones_curr_state =sum(MatStateID(i,:),2)*d;
     P_curr =P(i,:);
     check_2 = cost_prev.*(P_curr); % getting the expected value based on the new cost
     s_cost = check_2 *sun_ones_curr_state; % here i am not sure if we need to add the d or to multiply by the d
     ss_cost =sum(s_cost);
    cost(i) = min([cost(i), sum(ss_cost)]);
     if(delta<=cost(i))
    policy(h_, i) =1;
    end
    cost_prev = cost ;
end
end


ssiu =1;
%%
for h_=1:hbar
    for i  = 1:TotalStateAmount
     sss = sum(P(i,:) *sum(MatStateID(i,:),2)*d);
    cost(i) = min([cost(i), sum(P(i,:) *sum(MatStateID(i,:),2)*d)]);
    if(delta<=cost(i))
    policy(h_, i) =1;
    end
    end
end
ss = 2;
