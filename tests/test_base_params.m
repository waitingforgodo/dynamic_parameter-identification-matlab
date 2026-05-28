function test_base_params()
% 验证：最小参数集（base parameters）和原全参数集（full parameters）计算出来的力矩完全一样
% 用全参数算一遍力矩 τ_full
% 用最小参数算一遍力矩 τ_base
% 如果两者几乎一样（误差 < 1e-6）
% → 你的最小参数求解 完全正确
% ---------------------------------------------------------------------
% Test if base parameters are found correctly by comparing
% torque prediction of the model with standard parameters with 
% the torque prediction of the model with base parameters
% ----------------------------------------------------------------------

path_to_urdf = 'Marvin.urdf';
ur10 = parse_urdf(path_to_urdf);
no_links = 7;

% Perform QR decompostions (use cache to speed up repeated runs)
include_motor_dynamics = 0;
cacheFile = 'baseQR_7dof.mat';
[~, baseQR] = base_params_qr(include_motor_dynamics, cacheFile);

bb = baseQR.numberOfBaseParameters;% 最小参数数量
E = baseQR.permutationMatrix; % 置换矩阵
beta = baseQR.beta;% 关联矩阵
includeMotorDynamics = baseQR.motorDynamicsIncluded;

if includeMotorDynamics
    no_link_params = 11;% 10个连杆惯性参数 + 1个电机惯量
    ur10.pi(end+1,:) = rand(1,no_links);
else
    no_link_params = 10; % 纯连杆：10个惯性参数/连杆
end
% ur10.pi = reshape(ur10.pi,[no_link_params*no_links, 1]);
ur10.pi =ur10.pi(:);
% 打印 ur10.pi 的维度和长度
fprintf('ur10.pi size: [%d, %d], total: %d\n', size(ur10.pi,1), size(ur10.pi,2), numel(ur10.pi));

% Position, velocity and acceleration limits
q_min = -pi*ones(7,1);
q_max = pi*ones(7,1);
qd_max = 3*pi*ones(7,1);
q2d_max = 6*pi*ones(7,1);

% On random positions, velocities, aceeleations
for i = 1:25
    q_rnd = q_min + (q_max - q_min).*rand(7,1);
    qd_rnd = -qd_max + 2*qd_max.*rand(7,1);
    q2d_rnd = -q2d_max + 2*q2d_max.*rand(7,1);
    
    if includeMotorDynamics
        Yi = regressorWithMotorDynamics(q_rnd,qd_rnd,q2d_rnd);
    else
        g=[0 0 9.8065]';
        Yi = standard_regressor_marvin(q_rnd, qd_rnd, q2d_rnd,g);
    end
    % 用【全参数】计算力矩（标准答案）
    tau_full = Yi*ur10.pi;
    % 用【最小参数】计算力矩（测试值）
    pi_lgr_base = [eye(bb) beta]*E'*ur10.pi;
    Y_base = Yi*E(:,1:bb);
    tau_base = Y_base*pi_lgr_base;
    assert(norm(tau_full - tau_base) < 1e-6);
end
fprintf("Rigid Body Base Dynamics Test - OK!\n");