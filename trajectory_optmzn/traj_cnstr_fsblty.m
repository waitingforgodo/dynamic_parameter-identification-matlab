function [c, ceq] = traj_cnstr_fsblty(opt_vars, traj_par, baseQR)
% --------------------------------------------------------------------
% 轨迹可行性 + 激励度 复合约束函数
% 功能：
%    1) 先调用基础约束，保证关节位置、速度、加速度不超限
%    2) 可选：增加“激励度约束” —— 要求信息矩阵对数行列式 ≥ 阈值
%    3) 用于可行性搜索/优化，保证轨迹既能跑、又能有效辨识参数
% 输入：
%    opt_vars   - 待优化的傅里叶系数向量
%    traj_par   - 轨迹参数结构体
%    baseQR     - 基参数QR分解结果（最小参数集）
% 输出：
%    c    - 不等式约束 c <= 0
%    ceq  - 等式约束（本函数无，返回空）
% --------------------------------------------------------------------

% ===================== 1. 先加载基础可行性约束 =====================
% 调用基础约束函数：关节位置、速度、加速度限位 + 平滑性等
% 这些是保证轨迹能实际运行的基本约束
[c_basic, ceq] = traj_cnstr(opt_vars, traj_par);

% 初始化约束数组 = 基础约束
c = c_basic;

% ===================== 2. 判断是否需要添加激励度约束 =====================
% 如果 traj_par 中没有定义 min_logdet（最小激励阈值）
% 则直接返回，不添加激励约束
if ~isfield(traj_par, 'min_logdet') || isempty(traj_par.min_logdet)
    return;
end

% ===================== 3. 从优化变量还原轨迹 =====================
opt_vars = opt_vars(:);
N = traj_par.N;                % 傅里叶级数阶数
wf = traj_par.wf;              % 基频
T = traj_par.T;                % 轨迹总时间

% 选择时间采样点（优先使用优化粗网格 t_opt，速度更快）
if isfield(traj_par, 't_opt') && ~isempty(traj_par.t_opt)
    t = traj_par.t_opt;
else
    t = traj_par.t;
end

t = reshape(t, 1, []);% 保证时间向量为行向量
n = size(traj_par.q0, 1);      % 关节自由度 n = 7
ab = reshape(opt_vars, [2*n, N]);  % 把优化变量 reshape 成系数矩阵
a = ab(1:n, :);                % 傅里叶系数 a（位置部分）
b = ab(n+1:2*n, :);            % 傅里叶系数 b（速度部分）

% 计算多项式系数（用于轨迹起始/结束平滑）
c_pol = getPolCoeffs(T, a, b, wf, N, traj_par.q0);

% 生成轨迹：角度 q、速度 qd、加速度 q2d
[q, qd, q2d] = mixed_traj(t, c_pol, a, b, wf, N);

% ===================== 4. 构建观测矩阵 W（回归矩阵） =====================
% E1：基参数映射矩阵（全参数 → 最小独立基参数）
E1 = baseQR.permutationMatrix(:, 1:baseQR.numberOfBaseParameters);

nFriction = 3 * n;
W = zeros(numel(t) * n, baseQR.numberOfBaseParameters + nFriction);

% 遍历所有时间点，逐点构建回归矩阵
for i = 1:numel(t)
    if baseQR.motorDynamicsIncluded
        Yi = regressorWithMotorDynamics(q(:, i), qd(:, i), q2d(:, i));
    else
        Yi = standard_regressor_marvin(q(:, i), qd(:, i), q2d(:, i));
    end

    % 构建当前时间块的观测矩阵：
    %    Yi * E1        → 基参数回归部分
    %    frictionRegressor → 摩擦参数回归部分
    W((i-1)*n + (1:n), :) = [Yi * E1, frictionRegressor(qd(:, i))];
end

% ===================== 5. 列归一化（避免数值量级差异） =====================
% 计算每列的 L2 范数（列向量模长）
col_norms = sqrt(sum(W.^2, 1));

% 避免除零，把非常小的范数设为 1
col_norms(col_norms < 1e-10) = 1;

% 归一化后的观测矩阵
Wn = W ./ col_norms;

% ===================== 6. 计算信息矩阵 F = W' * W + 小对角项 =====================
F = Wn' * Wn + 1e-9 * eye(size(Wn, 2));  % 小对角项保证矩阵正定

% ===================== 7. 计算 logdet(F) —— 轨迹激励度指标 =====================
% 乔列斯基分解，判断正定性
[R, p] = chol(F);
if p ~= 0
    % 分解失败 → 矩阵奇异 → 激励极差
    logdetF = -inf;
else
    logdetF = 2 * sum(log(diag(R) + eps));
end


% ===================== 8. 添加激励度不等式约束 =====================
% 约束：min_logdet - logdetF ≤ 0
% 等价于：logdetF ≥ min_logdet
c(end + 1, 1) = traj_par.min_logdet - logdetF;
end
