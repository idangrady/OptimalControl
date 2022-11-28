function [idx_v,idx_s] = idx_to_state_velo(idx, Mv)
% assume Vs is first
% so every num_vm, we increment vs by one
idx_v = floor(idx/ Mv);
idx_s  =  idx- idx_v*Mv;
idx_v  = idx_v+1;
end
