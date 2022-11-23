function [c] = timeoptimalpathspeed2(omegai,si,tau,deltav,L2,Na,Mv)
%parameters:
%omegai     vector,   given omegavector
%si...      vector,   given distance vector
%tau...     constant, stepsize for the time
%deltav...  constant, stepsize for speed
%L2...      constant, limit for velocity change
%Na...      constant
%Mv...      constant

%% first: set up grid
  %1st calculate limits
  L1 = deltav*Na/tau;
  %L2 is given
  L3 = Mv*deltav;

  %2nd set up grid
  %grip is in 3 dimensions, v-dimension, s-dimension, tau-dimension
  deltas = deltav * tau;
  sl = si(end); %the largest s element
  Ms = sl/deltas; %elements of s-grid
  %Mv = elements of v-grid is already given
  h = length(si);   %expected that this is possible

  %problem with Ms
  v = 0:deltav:Mv*deltav;
  s = 0:deltas:Ms*deltas;
  s_normed = 0:1:Ms;
  v_normed = 0:1:Mv;
  C = ones(h,(Mv+1)*(Ms+1),(Mv+1)*(Ms+1))*inf;
  M = zeros(h,(Mv+1)*(Ms+1),(Mv+1)*(Ms+1));

  %fill grid
  for k = 1:h
    for v_currentIdx = 1:Mv+1 
      for s_current_Idx = 1:Ms+1
        for v_nextIdx = 1:Mv+1 
          for s_nextIdx = 1:Ms+1
            s_normed_current = s_normed(s_current_Idx);
            s_normed_next    = s_normed(s_nextIdx);
            v_normed_current = v_normed(v_currentIdx);
            v_normed_next    = v_normed(v_nextIdx);
            c2 = ((s_normed_next-s_normed_current-v_normed_current)*deltas - (v_normed_next-v_normed_current)*deltav*tau/2)/(tau^3*(1/6-1/4));
            c1 = (v_normed_next-v_normed_current)*deltav/tau-c2*tau/2;

            if((c2<0) && (c1>0))
              vintmax = v_normed_current-c1^2/2/c2;
            else
              vintmax = v_normed_current;
            end
          
            if(abs(c1)>L1)
              cost = inf;
            elseif(abs(c1+c2*tau)>L1)
              cost = inf;
            elseif(omegai(k)*vintmax^2 > L2)
              cost = inf;
            elseif(vintmax>L3)
              cost = inf;
            else
              cost = (si(k) - s_normed_current*deltas)^2; %assumption about quadratic error
            end
   
            %beginning contraint
            if(k==1) %s(0)=0 && %v(0)=0
              if(not((s_normed_current==0) && (v_normed_current==0)))
                cost = inf;
              end
            end
            
            %end contraint
            if(k==h) %s(h)=sl %v(0)=0
              if(not(v_normed_current==0))
                cost = inf;
              end
            end
          
            %cannot go backwards
            if(s_normed_next<s_normed_current)
              cost = inf;
            end
           
            %assign cost to grid
            C(k,s_current_Idx+Ms*(v_currentIdx-1),s_nextIdx+Ms*(v_nextIdx-1)) = cost;
            M(k,s_current_Idx+Ms*(v_currentIdx-1),s_nextIdx+Ms*(v_nextIdx-1)) = s_nextIdx+Ms*(v_nextIdx-1);          
          end
        end
      end
    end
  end

%% second: calculate optimal paths from grid
  % terminal cost
  number_of_elements = Mv*Ms;
  T = zeros(1,number_of_elements);
  size_M = size(M);
  h = size_M(1);
  for l = 1:length(T)
    J_{h+1}{l} = T(l);
  end
  for k = h:-1:1
    ni = number_of_elements; % state dimension
    for i=1:ni
      nj = number_of_elements;
      caux = zeros(1,nj);
      for j = 1:nj
        caux(j) = C(k,i,j) + J_{k+1}{M(k,i,j)};
      end
      [a,b] = sort(caux);
      J_{k}{i} = a(1); J{k}{i} = J_{k}{i};
      u{k}{i}(1) = b(1);
      for ell = 2:length(a)
        if abs( a(ell) - a(1))<1e-8
          u{k}{i}(ell) = b(ell);
        else
        break;
        end
      end
    end
  end

%% third, extract one optimal path from solutions
  len = length(J);
  u_one = zeros(len,1);
  J_one = zeros(len,1);
  
  J_tmp = cell2mat(J{1});
  J_one(1) = J_tmp(1);
  u_tmp = cell2mat(u{1});
  u_one(1) = u_tmp(1);
  
  for i =2:len
    J_tmp = cell2mat(J{i});
    u_tmp = cell2mat(u{i});
    J_one(i) = J_tmp(u_one(i-1));
    u_one(i) = u_tmp(u_one(i-1));
  end

%% fourth: extraxt optimal path and optimal speed
  len = length(u_one);
  s_one = zeros(len,1);
  v_one = zeros(len,1);
  
  s_one(1) = 0;
  v_one(1) = 0;
  
  for i = 2:len
    idx_v = floor((u_one(i-1)-1)/Ms)+1;
    idx_s = u_one(i-1)-Ms*(idx_v-1);
    v_one(i) = v(idx_v);
    s_one(i) = s(idx_s);
  end
  
 %% fifth: calculated optimal c from optimal speed and path

  c = zeros(2,length(v_one));
  for i = 1:length(v_one)
    if(i<length(v_one))
      v_normed_current = v_one(i)/deltav;
      v_normed_next = v_one(i+1)/deltav;
      s_normed_current = s_one(i)/deltas;
      s_normed_next = s_one(i+1)/deltas;
      c(2,i) = ((s_normed_next-s_normed_current-v_normed_current)*deltas - (v_normed_next-v_normed_current)*deltav*tau/2)/(tau^3*(1/6-1/4));
      c(1,i) =(v_normed_next-v_normed_current)*deltav/tau-c(2,i)*tau/2;
    else
      c(1,i) = c(1,i-1);
      c(2,i) = c(2,i-1);
    end
  end

figure()
hold on
plot(v_one,"g")
plot(s_one,"b-")
plot(si,"rx")
grid on;
hold off
legend("speed","optimal distance","wished distance")
  
end

