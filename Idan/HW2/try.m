

r_11=0.6; r_00 =0.2;r_01 = 0.2;
r={[1-r_00, 1-r_01; 1-r_01, 1-r_11], [r_00, r_01; r_01, r_11]};


id=0:2^N-1; %index numbers
n=size(r{1}, 1); %size of r
id_bin=dec2bin(id); % Convert index number in decimal integer to its binary representation


