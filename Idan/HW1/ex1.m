% Exercise notes

% state space Ms*r , Mv*r

% for each state si, vi we assign a cost function. -> make sure you
% determine the cost wisely such that the minimization of the function,
% would result in the min of the problem.

% loop over the values, the set the values which violate the constraints to
% inf

% parameters which will be sent via funciton 

Ms = 10;
Mv = 10;
r = 5;
delta_V = 1;
Na =1; %given param

delta_S = delta_V* r; %delta s

% create grid
Grid = zeros(Ms*r, Mv*r); % empty array for state space (Ms*r,Mv*r)

% create V state, S state
V= linspace(0, Mv, Mv+1) ;
S= linspace(0, Ms, Ms+1);

% create velues based on given params
L1 = (Na*delta_V)/r;
L3 = delta_V*Mv;
S_L = delta_S*Ms;

for vIdx = 1:Mv*r % the range is probably a mistake. I think it should not go from 1: Mv*r but perhaps via max steps?
    for sIdx = 1:Ms*r
        

        C_1_k = (delta_S*(6*(S(sIdx+1)- S(sIdx)- V(vIdx)))/(r^2) - 2*delta_V/r *(V(vIdx+1)-V(vIdx)));
        C_2_k = 3*delta_V/r^3 *(V(vIdx+1)-V(vIdx)) - 12*delta_S/r^3 *(V(sIdx+1)-S(sIdx)- V(vIdx));
        
        % 
        vint_max = V(vIdx); % understand whether the this K

         if (C_2_k<0)
             if  C_1_k>0
                 % understand what exactly is the k!
                 vint_max = vint_max- (C_1_k^2)/(2*C_2_k);
             end 
         end

        
        % constraint -> to check before assignment - to speed up process.

         if abs(C_1_k + C_2_k*r) >L1
            Grid(vIdx, sIdx)=inf;
        elseif abs(C_1_k) >L1
            Grid(vIdx, sIdx)=inf;

            % the last if statement is incorrect, need to add the 
         elseif (vint_max>L3)
           %constraint 1: V_intmax < delta_v *Mv else: inf
            Grid(vIdx, sIdx)=inf;
         
            %TODO:  add          
            %constraint : a(t) =c1_k+c2_k(t-k*r) < L2 (L2 is given) ; else inf
         
         else
         % cost function c1k+ c2k
            cost = C_1_k + C_2_k ; 
            Grid(vIdx, sIdx) = cost;
        end
        
       
    end
end