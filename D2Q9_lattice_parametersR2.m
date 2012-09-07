% 

function [w, ex, ey,bb_spd,Q] = D2Q9_lattice_parametersR2()

w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
ex = [0, 1, 0, -1, 0, 1, -1, -1, 1];
ey = [0, 0, 1, 0, -1, 1, 1, -1, -1];
bb_spd = [1, 4, 5, 2, 3, 8, 9, 6, 7];

Q = zeros(2*2,9);

c_s_sq = 1/3;

for spd=1:9
   c = [ex(spd) ey(spd)];
   Qtmp = c'*c - c_s_sq*eye(2);
   Q(:,spd)=Qtmp(:);
   % each row of Q represents the components of Q (from Latt dissertation)
  % Q(:,:,spd)= c'*c - c_s_sq*eye(2);
end