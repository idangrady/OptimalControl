
Mv  = 7;
Ms = 7;


idx = 1;

%[r,c] = ind2sub(Mv,idx) idx_s = mod(idx,Mv);


idx_v = floor(idx/ Mv);
idx_s  =  idx- idx_v*Mv;%max(idx - Mv*(idx_v-1),0);
idx_v  = idx_v+1;
s = [idx_v,idx_s]