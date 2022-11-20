function [cost_grid] = define_grid(omegai,si,tau,deltav,L2,Na,Mv)
%parameters:
%omegai     vector,   given omegavector
%si...      vector,   given distance vector
%tau...     constant, stepsize for the time
%deltav...  constant, stepsize for speed
%L2...      constant, limit for velocity change
%Na...      constant
%Mv...      constant


  %1st calculate limits
  L1 = deltav*Na/tau;
  %L2 is given
  L3 = Mv*deltav;



  %2nd set up grid
  %grip is in 3 dimensions, v-dimension, s-dimension, tau-dimension

  deltas = deltav * tau;
  sl = si(end); %the largest s element
  Ms = sl/deltas +1; %elements of s-grid
  %Mv = elements of v-grid is already given
  h = length(si);   %todo, assumption made


  cost_grid = zeros((Mv*Ms*Mv*Ms),h);    %Mv*Ms possible combinations of s & v, h is number of steps
  % v0,s0;
  % v0,s1
  %..
  %v0,sl
  %v1,s0
  %...
  %vmax,sl

  %problem with Ms
  v = linspace(0, Mv*deltav, Mv) ;
  s = linspace(0, Ms*deltas, Ms);

  %fill grid
  for k = 1:h
    for v_currentIdx = 1:Mv 
      for s_current_Idx = 1:Ms
        for v_nextIdx = 1:Mv 
          for s_nextIdx = 1:Ms

            s_current = s(s_current_Idx);
            s_next    = s(s_nextIdx);
            v_current = v(v_currentIdx);
            v_next    = v(v_nextIdx);

            %todo - think about initial state

            c2 = ((s_next-s_current-v_current)*deltas - (v_next-v_current)*deltav*tau/2)/(tau^3*(1/6-1/4));
            c1 = (v_next-v_current)*deltav/tau-c2*tau/2;

            if((c2<0) && (c1>0))
              vintmax = v_current-c1^2/2/c2;
            else
              vintmax = v_current;
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
  %             cost = c1 + c2;   %maybe wrong way to define the cost
              cost = (si(k) - s_current)^2; %assumption about quadratic error
            end
            
            %beginning contraint
            if(k==1) %s(0)=0 && %v(0)=0
              if(not((s_current==0) && (v_current==0)))
                cost = inf;
              end
            end
            
            %end contraint
            if(k==h) %s(h)=sl %v(0)=0
              if(not((s_current==sl) && (v_current==0)))
                cost = inf;
              end
            end
            %assign cost to grid
            cost_grid(s_nextIdx+(v_nextIdx-1)*Mv+(s_current_Idx-1)*(Ms*Mv)+ (v_currentIdx-1) * (Ms*Mv*Ms),k) = cost;
          end
        end
      end
    end
  end
  % v[k+1] = v[k] + 1/deltav *(c[1,k]*tau + c[2,k] *tau*tau/2);
  % s[k+1] = s[k] + v[k] + 1/deltas *(c[1,k]*tau*tau/2 + c[2,k] *tau*tau*tau/6);

  %       c2[k] = (s(k+1)-s[k]-v[k])*deltas - (v[k+1]-v[k])*deltav*tau/2)/(tau*tau*tau*(1/6-1/4));
  %       c1[k] = (v[k+1]-v[k])*deltav/tau-c2[k]*tau/2;

%include contrains of physical model


end

