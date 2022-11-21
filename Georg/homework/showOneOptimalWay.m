function [u_one,J_one] = showOneOptimalWay(u,J)
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
  
end

