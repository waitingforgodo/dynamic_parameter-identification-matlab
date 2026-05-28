function [c, ceq] = traj_cnstr(opt_vars, traj_par)
% -------------------------------------------------------------------------
% 轨迹优化：关节位置/速度/加速度 约束函数（归一化版本）
% 作用：给优化器提供不等式约束，保证生成的轨迹安全、可执行
% 工程优化点：
%    1) 支持使用密网格检查约束，保证全程不超限
%    2) 所有约束归一化到 [-1,1] 量级，让优化器更稳定、收敛更快
% 输入：
%    opt_vars ：优化变量（傅里叶系数 a + b 展平成的向量）
%    traj_par ：轨迹参数结构体（包含限位、时间网格、轨迹配置等）
% 输出：
%    c   ：不等式约束 c ≤ 0（优化器要求）
%    ceq ：等式约束（此处无，为空）
% -------------------------------------------------------------------------

%% 1. 读取轨迹基本参数
opt_vars = opt_vars(:);
N = traj_par.N;                % 傅里叶级数的谐波阶数
wf = traj_par.wf;              % 轨迹基频
T = traj_par.T;                % 轨迹总时长

%% 2. 选择约束检查的时间网格（优先使用密网格，保证安全）
if isfield(traj_par, 't_cnstr') && ~isempty(traj_par.t_cnstr)
    t = traj_par.t_cnstr;
else
    t = traj_par.t;
end

t = reshape(t, 1, []);
%% 3. 从优化变量中恢复轨迹系数
n = size(traj_par.q0, 1);                     % 机器人自由度（7）
ab = reshape(opt_vars, [2*n, N]);              % 把一维优化变量 reshape 成系数矩阵
a = ab(1:n, :);                               % 傅里叶系数 a（位置项）
b = ab(n+1:2*n, :);                           % 傅里叶系数 b（速度项）

%% 4. 计算轨迹多项式系数 + 生成完整轨迹
c_pol = getPolCoeffs(T, a, b, wf, N, traj_par.q0);% 计算多项式系数（保证起止点平滑）
[q, qd, q2d] = mixed_traj(t, c_pol, a, b, wf, N);% 生成轨迹：位置、速度、加速度

%% 5. 计算归一化分母（防止除零，统一约束量级）
% 位置范围（最大-最小），用于归一化位置约束
q_span = max(traj_par.q_max(:) - traj_par.q_min(:), 1e-6);
% 速度/加速度上限，用于归一化速度/加速度约束
qd_lim = max(traj_par.qd_max(:), 1e-6);
q2d_lim = max(traj_par.q2d_max(:), 1e-6);

%% 6. 构建 4 组归一化不等式约束（全部要求 c ≤ 0）
c = zeros(4*n, 1);  % 4类约束 × n个关节 = 4n个约束值

% 约束1：各关节【最小位置】 ≥ q_min
% 公式：(q_min - min(q)) ≤ 0 → 归一化
c(1:n)       = (traj_par.q_min(:) - min(q, [], 2)) ./ q_span;

% 约束2：各关节【最大位置】 ≤ q_max
% 公式：(max(q) - q_max) ≤ 0 → 归一化
c(n+1:2*n)   = (max(q, [], 2) - traj_par.q_max(:)) ./ q_span;

% 约束3：各关节【最大绝对速度】 ≤ qd_max
% 公式：(max|qd| - qd_max) ≤ 0 → 归一化
c(2*n+1:3*n) = (max(abs(qd), [], 2) - traj_par.qd_max(:)) ./ qd_lim;

% 约束4：各关节【最大绝对加速度】 ≤ q2d_max
% 公式：(max|q2d| - q2d_max) ≤ 0 → 归一化
c(3*n+1:4*n) = (max(abs(q2d), [], 2) - traj_par.q2d_max(:)) ./ q2d_lim;

%% 7. 无等式约束
ceq = [];
end
