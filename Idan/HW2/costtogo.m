
P =[1/3 1/12 1/6 1/12 1/3];
Val = [8 9 10 11 12];

cost = Val;
cost_prev = Val;

for h_ =1:2
    for i =1:length(cost)
    sum_cost = sum(cost_prev.*P);
    cost(i) = min([cost(i), sum(sum_cost)]);
    end
    cost_prev = cost ;
end

s = 2;