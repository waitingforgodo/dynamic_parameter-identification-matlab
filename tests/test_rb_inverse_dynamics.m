function test_rb_inverse_dynamics()
% 验证你自己写的机器人动力学模型（质量矩阵 M、科氏力 C、重力 G、回归矩阵 Y）完全正确，和 MATLAB 官方工具箱结果一致。
path_to_urdf = 'Marvin.urdf';
ur10 = parse_urdf(path_to_urdf);

rbt = importrobot('Marvin.urdf');
rbt.DataFormat = 'column';
rbt.Gravity = [0 0 -9.81];

no_iter = 100;
for i = 1:no_iter
    q = -2*pi + 4*pi*rand(7,1);
    q_d = zeros(7,1);
    q_2d = zeros(7,1);

    Ylgr = standard_regressor_marvin(q,q_d,q_2d);

    tau_matlab = inverseDynamics(rbt,q,q_d,q_2d);
    tau_reg = Ylgr*reshape(ur10.pi,[70,1]);
    tau_manip = M_mtrx_fcn(q, ur10.pi(:))*q_2d + ...
                C_mtrx_fcn(q, q_d, ur10.pi(:))*q_d + ...
                G_vctr_fcn(q, ur10.pi(:));

%   verifying if regressor is computed correctly
    assert(norm(tau_matlab - tau_reg) < 1e-8);
    assert(norm(tau_matlab - tau_manip) < 1e-8);
end

fprintf("Rigid Body Inverse Dynamics Test - OK!\n");