
% prepares a vector representing destination matrix entries during a fully
% periodic streaming operation.

% simple, but slow...and that's just fine.

function target_vec = genTargetVecD2Q9r2(lx,ly)


target_mat = zeros(lx*ly,9);
nnode = lx*ly;

for x = 1:lx
    for y = 1:ly
        yn = mod(y,ly)+1; 
        xe = mod(x,lx)+1;
        ys = mod(y-2,ly)+1;
        xw = mod(x-2,lx)+1;
        source_dof = lx*(y-1)+x;
        
        % stationary particle
        dir = 1;
        dof_num = (lx*(y-1)+x);
        f_num = nnode*(dir-1)+dof_num;
        target_mat(source_dof,dir) = f_num;
        
        % east particle
        dir = 2;
        dof_num = lx*(y-1)+xe;
        f_num = nnode*(dir-1)+dof_num;
        target_mat(source_dof,dir) = f_num;
        
        %north particle
        dir = 3;
        dof_num = lx*(yn-1) + x;
        f_num = nnode*(dir-1) + dof_num;
        target_mat(source_dof,dir) = f_num;
        
        %west particle
        dir = 4;
        dof_num = lx*(y-1)+xw;
        f_num = nnode*(dir-1) + dof_num;
        target_mat(source_dof,dir) = f_num;
        
        %south particle
        dir = 5;
        dof_num = lx*(ys-1)+x;
        f_num = nnode*(dir-1)+dof_num;
        target_mat(source_dof,dir) = f_num;
        
        %north-east particle
        dir = 6;
        dof_num = lx*(yn-1)+xe;
        f_num = nnode*(dir-1)+dof_num;
        target_mat(source_dof,dir) = f_num;
        
        %north-west particle
        dir = 7;
        dof_num = lx*(yn - 1) + xw;
        f_num = nnode*(dir-1)+dof_num;
        target_mat(source_dof,dir) = f_num;
        
        %south-west particle
        dir = 8;
        dof_num = lx*(ys - 1)+xw;
        f_num = nnode*(dir-1)+dof_num;
        target_mat(source_dof,dir) = f_num;
        
        %south-east particle
        dir = 9;
        dof_num = lx*(ys - 1) + xe;
        f_num = nnode*(dir-1)+dof_num;
        target_mat(source_dof,dir) = f_num;
        
    end
end

target_vec = target_mat(:);





