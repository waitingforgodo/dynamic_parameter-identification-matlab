function out = traj_cost_lgr(opt_vars,traj_par,baseQR,model_type)
% -------------------------------------------------------------------
% 轨迹优化代价函数（进阶工程版）
% 目标：评估一条轨迹是否适合动力学参数辨识（越小越好）
% -------------------------------------------------------------------

%% 1. 初始化与变量解析
% ---------------------------------------------------------------
opt_vars = opt_vars(:);
n  = size(traj_par.q0, 1); % 关节数量

N  = traj_par.N;
wf = traj_par.wf;
T  = traj_par.T;

if isfield(traj_par, 't_opt') && ~isempty(traj_par.t_opt)
    t = traj_par.t_opt;
else
    t = traj_par.t;
end
t = reshape(t, 1, []);
    
% 解析傅里叶系数 (a: sin, b: cos)
ab = reshape(opt_vars, [2*n, N]);
a = ab(1:n, :);
b = ab(n+1:2*n, :);

%% 2. 轨迹生成与物理约束检查
% ---------------------------------------------------------------
% 计算多项式系数并生成混合轨迹 (位置 q, 速度 qd, 加速度 q2d)
c_pol = getPolCoeffs(T, a, b, wf, N, traj_par.q0);
[q, qd, q2d] = mixed_traj(t, c_pol, a, b, wf, N);

% --- 检查硬性约束（越限惩罚） ---
viol = normalizedLimitViolation(q, qd, q2d, traj_par);
viol_pos = max(viol, 0); % 只保留正值违规量

% 如果存在任何违规，返回巨大惩罚值，迫使优化器远离该区域
% if any(viol_pos > 0)
    % % 计算运动强度作为辅助参考
    % motion_level = mean(std(qd, 0, 2) ./ max(traj_par.qd_max(:), 1e-6));
    % out = 1e6 * (1 + sum(viol_pos.^2)) + 1e3 / (motion_level + 1e-6);
    % return;
% end

% 惩罚系数设为 1e4，既能起到限制作用，又不会掩盖激励指标本身的数值变化
constraint_penalty = 1e4 * sum(viol_pos.^2); 


%% 3. 构建观测矩阵 (回归矩阵 W)
% ---------------------------------------------------------------
E1 = baseQR.permutationMatrix(:, 1:baseQR.numberOfBaseParameters);
nFriction = 3 * n; 

% 预分配内存，提高循环效率
W = zeros(numel(t) * n, baseQR.numberOfBaseParameters + nFriction);

for i = 1:numel(t)
    if strcmpi(model_type, 'Newton-Euler')
        Yi = SymForm_NT_Y(q(:,i),qd(:,i),q2d(:,i));
    elseif strcmpi(model_type, 'Lagrangian')
        Yi = SymForm_LG_Y(q(:,i),qd(:,i),q2d(:,i));
    else
        assert(false,'Model type error');
    end
    
    Fi = SymForm_F(qd(:, i));

    % 组装当前时刻的回归矩阵行
    row_idx = (i-1)*n + (1:n);
    W(row_idx, :) = [Yi * E1, Fi];
end

% 异常值检查
if any(~isfinite(W(:)))
    out = 1e9;
    return;
end

%% 4. 计算激励指标 (基于归一化信息矩阵)
% ---------------------------------------------------------------
% 对 W 的每一列进行归一化，消除不同参数量级差异的影响
col_norms = sqrt(sum(W.^2, 1));
col_norms(col_norms < 1e-10) = 1; % 防止除以零
Wn = W ./ col_norms;

% 计算 Fisher 信息矩阵 F
F = Wn' * Wn + 1e-9 * eye(size(Wn, 2));

if any(~isfinite(F(:)))
    out = 1e9;
    return;
end

% 使用 Cholesky 分解计算 logdet(F)，数值稳定性更好
[R, p] = chol(F);
if p ~= 0
    % 如果矩阵不正定，给予较大惩罚
    out = 1e8 + trace(F);
    return;
end
logdetF = 2 * sum(log(diag(R) + eps));


%% 5. 综合代价函数计算
% ---------------------------------------------------------------
% A. 运动不足惩罚：鼓励轨迹动起来，避免静止或微动
motion_level = mean(std(qd, 0, 2) ./ max(traj_par.qd_max(:), 1e-6));
motion_penalty = 2.0 / (motion_level + 1e-3);

% B. 正则化惩罚：限制优化变量幅值，防止系数过大导致震荡
regularization = 1e-5 * (opt_vars' * opt_vars) / numel(opt_vars);

% C. 对称性惩罚：鼓励正负方向运动幅度均衡，全面激发摩擦参数
q_above = max(q - traj_par.q0, [], 2);
q_below = max(traj_par.q0 - q, [], 2);
symmetry_ratio = min(q_above, q_below) ./ (max(q_above, q_below) + 1e-6);
symmetry_penalty = 55.0 * sum(1 - symmetry_ratio);

% --- 最终输出：代价 = -激励强度 + 各类惩罚项 ---
% out = -logdetF + motion_penalty + symmetry_penalty + regularization;
out = -logdetF + motion_penalty + symmetry_penalty + regularization + constraint_penalty;

end


%% 辅助子函数：计算归一化约束违反程度
function viol = normalizedLimitViolation(q, qd, q2d, traj_par)
    q_span = max(traj_par.q_max(:) - traj_par.q_min(:), 1e-6);
    qd_lim = max(traj_par.qd_max(:), 1e-6);
    q2d_lim = max(traj_par.q2d_max(:), 1e-6);
    
    % 分别计算位置下限、上限，以及速度和加速度的越限情况
    viol = [...
        (traj_par.q_min(:) - min(q, [], 2)) ./ q_span; ... % 位置下限违规
        (max(q, [], 2) - traj_par.q_max(:)) ./ q_span; ... % 位置上限违规
        (max(abs(qd), [], 2) - traj_par.qd_max(:)) ./ qd_lim; ... % 速度违规
        (max(abs(q2d), [], 2) - traj_par.q2d_max(:)) ./ q2d_lim]; % 加速度违规
end




