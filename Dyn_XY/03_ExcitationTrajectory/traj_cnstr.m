function [c,ceq] = traj_cnstr(opt_vars,traj_par,baseQR,model_type)
% -------------------------------------------------------------------------
% 轨迹优化复合约束函数（可行性 + 激励度）
% 功能：
%    1. 基础约束：保证关节位置、速度、加速度不超限（归一化处理）
%    2. 进阶约束：当传入 baseQR 且定义了 min_logdet 时，增加激励度约束
% 输入：
%    opt_vars ：优化变量（傅里叶系数 a + b 展平成的向量）
%    traj_par ：轨迹参数结构体（包含限位、时间网格、轨迹配置等）
%    baseQR   ：(可选) 基参数QR分解结果结构体。若为空或不传，则只进行可行性约束
% 输出：
%    c   ：不等式约束 c ≤ 0
%    ceq ：等式约束（此处无，为空）
% -------------------------------------------------------------------------


%% 1. 初始化与参数解析
opt_vars = opt_vars(:);
n = size(traj_par.q0, 1); % 机器人自由度 

N = traj_par.N;     % 傅里叶级数谐波阶数
wf = traj_par.wf;   % 轨迹基频
T = traj_par.T;     % 轨迹总时长

% 恢复傅里叶系数
ab = reshape(opt_vars, [2*n, N]);
a = ab(1:n, :);                               % 位置项系数
b = ab(n+1:2*n, :);                           % 速度项系数

%% 2. 物理限位约束 (使用密网格 t_cnstr 保证安全)
if isfield(traj_par, 't_cnstr') && ~isempty(traj_par.t_cnstr)
    t_cnstr = traj_par.t_cnstr;
else
    t_cnstr = traj_par.t;
end
t_cnstr = reshape(t_cnstr, 1, []);

% 生成用于检查物理限位的轨迹
c_pol = getPolCoeffs(T, a, b, wf, N, traj_par.q0);
[q, qd, q2d] = mixed_traj(t_cnstr, c_pol, a, b, wf, N);

% 计算归一化分母（防止除零，统一量级）
q_span = max(traj_par.q_max(:) - traj_par.q_min(:), 1e-6);
qd_lim = max(traj_par.qd_max(:), 1e-6);
q2d_lim = max(traj_par.q2d_max(:), 1e-6);

% 构建前 4n 个基础不等式约束
c = zeros(4*n, 1);
c(1:n)       = (traj_par.q_min(:) - min(q, [], 2)) ./ q_span;               % 最小位置
c(n+1:2*n)   = (max(q, [], 2) - traj_par.q_max(:)) ./ q_span;               % 最大位置
c(2*n+1:3*n) = (max(abs(qd), [], 2) - traj_par.qd_max(:)) ./ qd_lim;        % 最大速度
c(3*n+1:4*n) = (max(abs(q2d), [], 2) - traj_par.q2d_max(:)) ./ q2d_lim;     % 最大加速度

%% 3. 激励度约束 (仅在提供 baseQR 和 min_logdet 时激活)
if nargin >= 3 && ~isempty(baseQR) && isfield(traj_par, 'min_logdet') && ~isempty(traj_par.min_logdet)
    
    % 优先使用专门用于优化的粗网格 t_opt 以提升计算速度
    if isfield(traj_par, 't_opt') && ~isempty(traj_par.t_opt)
        t_opt = traj_par.t_opt;
    else
        t_opt = traj_par.t;
    end
    t_opt = reshape(t_opt, 1, []);
    
    % 重新生成用于构建 W 矩阵的轨迹 (如果 t_opt 和 t_cnstr 相同，这里可进一步优化缓存)
    [q_opt, qd_opt, q2d_opt] = mixed_traj(t_opt, c_pol, a, b, wf, N);
    
    % 提取基参数映射矩阵
    E1 = baseQR.permutationMatrix(:, 1:baseQR.numberOfBaseParameters);
    nFriction = 2 * n; % 库伦+粘滞摩擦模型，每个关节2个参数
    nTime = numel(t_opt);
    
    W = zeros(nTime * n, baseQR.numberOfBaseParameters + nFriction);
    
    % 逐点构建观测矩阵 W
    for i = 1:nTime
        if strcmpi(model_type, 'Newton-Euler')
            Yi = SymForm_NT_Y(q_opt(:,i),qd_opt(:,i),q2d_opt(:,i));
        elseif strcmpi(model_type, 'Lagrangian')
            Yi = SymForm_LG_Y(q_opt(:,i),qd_opt(:,i),q2d_opt(:,i));
        else
            assert(false,'Model type error');
        end
        % 拼接动力学回归项与摩擦回归项
        W((i-1)*n + (1:n), :) = [Yi * E1, SymForm_F(qd_opt(:,i))];
    end
    
    % 列归一化 (消除不同物理量纲的影响)
    col_norms = sqrt(sum(W.^2, 1));
    col_norms(col_norms < 1e-10) = 1;
    Wn = W ./ col_norms;
    
    % 计算信息矩阵 F 及其对数行列式 logdet(F)
    F = Wn' * Wn + 1e-9 * eye(size(Wn, 2));
    [R, p] = chol(F);
    if p ~= 0
        logdetF = -inf; % 矩阵奇异，给予极差的惩罚值
    else
        logdetF = 2 * sum(log(diag(R) + eps));
    end
    
    % 添加激励度约束：要求 logdetF >= min_logdet
    c(end + 1, 1) = traj_par.min_logdet - logdetF;
end

%% 4. 无等式约束
ceq = [];
end