function Y = regressorWithMotorDynamics(q,qd,q2d)
% ----------------------------------------------------------------------
% This function adds motor dynamics to rigid body regressor.
% It is simplified model of motor dynamics, it adds only reflected
% inertia i.e. I_rflctd = Im*N^2 where N is reduction ratio - I_rflctd*q_2d
% parameter is added to existing vector of each link [pi_i I_rflctd_i]
% so that each 7-DOF link has 11 parameters (total 77 parameters)
% ----------------------------------------------------------------------
% 输入校验：修改为 7自由度 7x1 向量
if size(q,1)==7 && size(q,2)==1 && size(qd,1)==7 && size(qd,2)==1 ...
        && size(q2d,1)==7 && size(q2d,2)==1
    % 调用 7自由度 刚体回归器（请替换为你自己的7自由度函数名）
    Y_rgd_bdy = standard_regressor_marvin(q,qd,q2d); 
    
    % 电机动力学：7x7 对角矩阵（角加速度）
    Y_mtrs = diag(q2d);
    
    % 7轴拼接：每10个刚体参数 + 1个电机惯量参数，共7组
    Y = [Y_rgd_bdy(:,1:10),  Y_mtrs(:,1), ...
         Y_rgd_bdy(:,11:20), Y_mtrs(:,2), ...
         Y_rgd_bdy(:,21:30), Y_mtrs(:,3), ...
         Y_rgd_bdy(:,31:40), Y_mtrs(:,4), ...
         Y_rgd_bdy(:,41:50), Y_mtrs(:,5), ...
         Y_rgd_bdy(:,51:60), Y_mtrs(:,6), ...
         Y_rgd_bdy(:,61:70), Y_mtrs(:,7)];
else
    error('Input dimension error! Please input 7x1 joint vectors.');
end