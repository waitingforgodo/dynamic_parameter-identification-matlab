% ---------------------------------------------------------------------
% 7自由度机械臂动力学参数辨识轨迹优化
% 工程优化版本特点：
% 1) 自动解析URDF路径，鲁棒性更强
% 2) 基于物理意义设置轨迹系数边界
% 3) 目标函数用粗网格加速，约束/验证用密网格保证精度
% 4) 自动搜索可行初始点，避免优化直接失败
% 5) 使用logdet激励指标作为优化目标，比直接优化条件数更稳定
% ---------------------------------------------------------------------
robot = parse_urdf('Marvin.urdf'); %#ok<NASGU>

include_motor_dynamics = 0;
% 基参数QR分解缓存文件，避免重复计算
cacheFile = fullfile('trajectory_optmzn', 'baseQR_7dof.mat');
% 计算基参数（最小参数集），返回基参数相关矩阵
[~, baseQR] = base_params_qr(include_motor_dynamics, cacheFile);

% 选择优化算法：patternsearch(模式搜索) / ga(遗传算法)
% optmznAlgorithm = 'patternsearch';
optmznAlgorithm = 'ga';

% Trajectory parameters
traj_par.n_dof = 7;
traj_par.T =60; % 轨迹总时间
traj_par.wf = 2 * pi / traj_par.T; % 基频率
traj_par.N =5; % 傅里叶级数阶数
traj_par.q0 = deg2rad([0 0 0 -42.5 0 0 0]');
traj_par.q_min = deg2rad([-170  -110  -170  -130  -170  -40  -40]');
traj_par.q_max = deg2rad([ 170   110   170    50   170   40   40]');
traj_par.qd_max = [3.0 3.0 3.0 3.0 3.0 3.0 3.0]';
traj_par.q2d_max = [110 100 110 110 110 110 110]';

% cond (W) = 平衡性（越小越好，难优化）
% logdetF = 激励强度（越大越好，易优化）
% logdetF 仅放入代价函数（-logdetF），不作为硬约束，否则很难找可行点
traj_par.min_logdet = [];
traj_par.bound_safety = 1.0;

%采样网格设置（工程优化：目标粗网格，约束密网格）
% 目标函数计算用粗网格，减少计算量，加速优化
% Separate grids: accelerate objective, keep constraints dense
traj_par.t_opt = 0:0.25:traj_par.T;
traj_par.t_cnstr = 0:0.1:traj_par.T;
% 绘图用更密网格
traj_par.t_plot = 0:0.05:traj_par.T;
traj_par.t = traj_par.t_cnstr;
% 时间采样步长
traj_par.t_smp = traj_par.t_cnstr(2) - traj_par.t_cnstr(1);
% 目标利用率 60%-70%（相对 qd_max / q2d_max 峰值）
traj_par.target_qd_util  = 0.65;
traj_par.target_q2d_util = 0.60;
traj_par.w_util_qd  = 20.0;
traj_par.w_util_q2d = 20.0;
traj_par.w_motion = 15.0;
traj_par.w_symmetry = 50.0;
% Validate q0
for i = 1:traj_par.n_dof
    if traj_par.q0(i) < traj_par.q_min(i) || traj_par.q0(i) > traj_par.q_max(i)
        warning('关节 %d 的初始位置 q0=%.2f° 超出限位 [%.2f°, %.2f°]', ...
            i, rad2deg(traj_par.q0(i)), rad2deg(traj_par.q_min(i)), rad2deg(traj_par.q_max(i)));
    end
end

% 待优化变量总个数：2*N*n_dof（a和b两组傅里叶系数）
n_dof = traj_par.n_dof;
n_opt = 2 * n_dof * traj_par.N;

% 线性不等式/等式约束（本程序无全局线性约束）
A = []; b = [];
Aeq = []; beq = [];

% 傅里叶系数搜索边界 lb/ub（位置限位由 traj_cnstr 非线性约束保证）
[lb, ub, traj_par.ab_bound] = build_fourier_bounds(traj_par, n_opt);
fprintf('Coeff bounds: per-joint max ~ [%.3f, %.3f, ...] rad/s\n', ...
    traj_par.ab_bound(1), traj_par.ab_bound(min(2, n_dof)));


%%自动搜索可行初始点（避免优化从不可行点开始直接失败）
[x0, seedInfo] = find_feasible_seed(lb, ub, traj_par, baseQR);
fprintf('Seed search: feasible=%d, best_violation=%.3e, best_cost=%.3f\n', ...
    seedInfo.feasibleFound, seedInfo.bestViolation, seedInfo.bestCost);

% 确保并行池已启动
if isempty(gcp('nocreate'))
    parpool;
end

if strcmp(optmznAlgorithm, 'patternsearch')
    optns_pttrnSrch = optimoptions('patternsearch');
    optns_pttrnSrch.Display = 'iter';
    optns_pttrnSrch.StepTolerance = 1e-3;% 步长终止阈值
    optns_pttrnSrch.MeshTolerance = 1e-3;% 网格终止阈值
    optns_pttrnSrch.FunctionTolerance = 1e-3;% 目标函数终止阈值
    optns_pttrnSrch.ConstraintTolerance = 1e-6;% 约束终止阈值
    optns_pttrnSrch.MaxTime = 1800;% 最大时间
    optns_pttrnSrch.MaxIterations = 250;% 最大迭代次数
    optns_pttrnSrch.MaxFunctionEvaluations = 8e4;% 最大函数评估次数
    optns_pttrnSrch.UseCompletePoll = true; % 关闭全方向搜索
    optns_pttrnSrch.UseCompleteSearch = false; % 关闭全局搜索
    optns_pttrnSrch.UseParallel=true;        % 开启并行
    optns_pttrnSrch.UseVectorized=false;        % 关闭向量化
    % 执行模式搜索优化
    % [x, fval, exitflag, output] = patternsearch(@(x) traj_cost_lgr(x, traj_par, baseQR), x0, ...
    %     A, b, Aeq, beq, lb, ub, @(x) traj_cnstr(x, traj_par), optns_pttrnSrch);

    % x        = 最优轨迹
    % fval     = 轨迹有多好（分数越低越好）
    % exitflag = 优化成功还是失败
    % output   = 迭代过程详情
    [x, fval, exitflag, output] = patternsearch(@(x) traj_cost_lgr(x, traj_par, baseQR), x0, ...
        A, b, Aeq, beq, lb, ub, @(x) traj_cnstr(x, traj_par), optns_pttrnSrch);

elseif strcmp(optmznAlgorithm, 'ga')
    optns_ga = optimoptions('ga');
    optns_ga.Display = 'iter';
    optns_ga.MaxGenerations = 40;% 最大迭代次数
    optns_ga.PopulationSize = 40;% 种群大小
    optns_ga.EliteCount = 8;% 精英个体数
    optns_ga.FunctionTolerance = 1e-4;
    optns_ga.InitialPopulationRange = [lb.'; ub.'];
    optns_ga.PlotFcn = {@gaplotbestf, @gaplotmaxconstr};
    optns_ga.UseParallel = true;
    optns_ga.UseVectorized = false;

    % 将 find_feasible_seed 找到的 x0 注入 GA 初始种群（第 1 个个体为种子）
    popSize = optns_ga.PopulationSize;
    initPop = repmat(lb', popSize, 1) + rand(popSize, n_opt) .* repmat((ub - lb)', popSize, 1);
    initPop(1, :) = min(max(x0', lb'), ub');
    optns_ga.InitialPopulationMatrix = initPop;

    % Hybrid patternsearch：GA结束后局部精调，配上严格终止条件防止跑不停
    hybridOpts = optimoptions('patternsearch');
    hybridOpts.MaxIterations = 100;
    hybridOpts.MaxFunctionEvaluations = 1e4;
    hybridOpts.MaxTime = 100;
    hybridOpts.StepTolerance = 1e-3;
    hybridOpts.FunctionTolerance = 1e-3;
    hybridOpts.UseCompletePoll = true;
    hybridOpts.UseParallel = true;
    hybridOpts.Display = 'iter';
    optns_ga.HybridFcn = {@patternsearch, hybridOpts};

    % 执行遗传算法优化
    % [x, fval, exitflag, output] = ga(@(x) traj_cost_lgr(x, traj_par, baseQR), n_opt, ...
    %     A, b, Aeq, beq, lb, ub, @(x) traj_cnstr(x, traj_par), optns_ga);

    % 执行遗传算法优化
    [x, fval, exitflag, output] = ga(@(x) traj_cost_lgr(x, traj_par, baseQR), n_opt, ...
        A, b, Aeq, beq, lb, ub, @(x) traj_cnstr(x, traj_par), optns_ga);
else
    error('Chosen algorithm is not found among implemented ones');
end

% [c_final, ~] = traj_cnstr(x, traj_par);
% fprintf('Stage-1 finished. Exitflag=%d, max normalized violation=%.3e\n', ...
%     exitflag, max(c_final));

[c_final, ~] = traj_cnstr(x, traj_par);
fprintf('Stage-1 finished. Exitflag=%d, max normalized violation=%.3e\n', ...
    exitflag, max(c_final));

% 可选：fmincon 局部精调（如果安装了优化工具箱）
refineInfo = struct('used', false, 'accepted', false, 'exitflag', [], 'output', []);
if exist('fmincon', 'file') == 2
    try
        optns_fmincon = optimoptions('fmincon');
        optns_fmincon.Display = 'iter';
        optns_fmincon.Algorithm = 'sqp'; % 序列二次规划
        optns_fmincon.MaxFunctionEvaluations = 2e4;
        optns_fmincon.MaxIterations = 120;% 最大迭代次数
        optns_fmincon.StepTolerance = 1e-4;% 步长终止阈值
        optns_fmincon.ConstraintTolerance = 1e-6;% 约束终止阈值
        optns_fmincon.OptimalityTolerance = 1e-4;

        refineInfo.used = true;
        % 局部精调优化
        % [x_ref, fval_ref, exitflag_ref, output_ref] = fmincon(@(x) traj_cost_lgr(x, traj_par, baseQR), x, ...
        %     A, b, Aeq, beq, lb, ub, @(x) traj_cnstr(x, traj_par), optns_fmincon);

        % 局部精调优化
        [x_ref, fval_ref, exitflag_ref, output_ref] = fmincon(@(x) traj_cost_lgr(x, traj_par, baseQR), x, ...
            A, b, Aeq, beq, lb, ub, @(x) traj_cnstr(x, traj_par), optns_fmincon);

        [c_ref, ~] = traj_cnstr(x_ref, traj_par);
        refineInfo.exitflag = exitflag_ref;
        refineInfo.output = output_ref;

        if all(c_ref <= 1e-6) && fval_ref < fval
            % 接受更优的精调结果
            x = x_ref;
            fval = fval_ref;
            exitflag = exitflag_ref;
            output = output_ref;
            c_final = c_ref;
            refineInfo.accepted = true;
        end
    catch ME
        warning('fmincon refinement skipped: %s', ME.message);
    end
end

fprintf('Optimization finished. Exitflag=%d, max normalized violation=%.3e, objective=%.6f\n', ...
    exitflag, max(c_final), fval);

% Plot optimized trajectory
ab = reshape(x, [2*n_dof, traj_par.N]);
a = ab(1:n_dof, :);
b = ab(n_dof+1:2*n_dof, :);
c_pol = getPolCoeffs(traj_par.T, a, b, traj_par.wf, traj_par.N, traj_par.q0);
[q, qd, q2d] = mixed_traj(traj_par.t_plot, c_pol, a, b, traj_par.wf, traj_par.N);

qd_peak = max(abs(qd), [], 2);
q2d_peak = max(abs(q2d), [], 2);
fprintf('\nPeak utilization (%% of limit):\n');
for j = 1:n_dof
    fprintf('  Joint %d: qd=%.1f%%, q2d=%.1f%%\n', j, ...
        100*qd_peak(j)/traj_par.qd_max(j), 100*q2d_peak(j)/traj_par.q2d_max(j));
end
fprintf('  Mean: qd=%.1f%%, q2d=%.1f%%\n\n', ...
    100*mean(qd_peak./traj_par.qd_max), 100*mean(q2d_peak./traj_par.q2d_max));

jointLabels = arrayfun(@(k) sprintf('q%d', k), 1:n_dof, 'UniformOutput', false);
qdLabels = arrayfun(@(k) sprintf('qd%d', k), 1:n_dof, 'UniformOutput', false);
q2dLabels = arrayfun(@(k) sprintf('q2d%d', k), 1:n_dof, 'UniformOutput', false);

figure
subplot(3,1,1)
plot(traj_par.t_plot, q)
ylabel('$q$', 'interpreter', 'latex')
grid on
legend(jointLabels{:})
subplot(3,1,2)
plot(traj_par.t_plot, qd)
ylabel('$\dot{q}$', 'interpreter', 'latex')
grid on
legend(qdLabels{:})
subplot(3,1,3)
plot(traj_par.t_plot, q2d)
ylabel('$\ddot{q}$', 'interpreter', 'latex')
grid on
legend(q2dLabels{:})

% Save optimized trajectory coefficients
pathToFolder = 'trajectory_optmzn/optimal_trjctrs/';
if ~exist(pathToFolder, 'dir')
    mkdir(pathToFolder);
end

t1 = strcat('N', num2str(traj_par.N), 'T', num2str(traj_par.T));
if strcmp(optmznAlgorithm, 'patternsearch')
    filename = strcat(pathToFolder, 'ptrnSrch_', t1, '_7dof_QR.mat');
elseif strcmp(optmznAlgorithm, 'ga')
    filename = strcat(pathToFolder, 'ga_', t1, '_7dof_QR.mat');
elseif strcmp(optmznAlgorithm, 'fmincon')
    filename = strcat(pathToFolder, 'fmncn_', t1, '_7dof_QR.mat');
end

save(filename, 'a', 'b', 'c_pol', 'traj_par', 'x', 'fval', 'exitflag', 'output', 'seedInfo', 'refineInfo');
fprintf('Saved trajectory to %s\n', filename);

%% ===================== 子函数：傅里叶系数搜索边界 =====================
function [lb, ub, ab_bound] = build_fourier_bounds(traj_par, n_opt)
% 为 GA/PS 提供 lb/ub 搜索盒。
% 仅依据速度、加速度估算系数上界；关节位置是否超限由 traj_cnstr 密网格约束检查。
% 若把位置也放进 min(...)，T=120 时低阶谐波会被 pos_margin/(k*wf) 严重压小。
N = traj_par.N;
n_dof = traj_par.n_dof;
wf = traj_par.wf;

if nargin < 2 || isempty(n_opt)
    n_opt = 2 * n_dof * N;
end

if isfield(traj_par, 'bound_safety') && ~isempty(traj_par.bound_safety)
    safety = traj_par.bound_safety;
else
    safety = 0.90;
end

lb = zeros(n_opt, 1);
ub = zeros(n_opt, 1);
bound_per_harmonic = zeros(n_dof, N);

for k = 1:N
    wk = k * wf;
    % 单谐波速度幅值上界：密网格 traj_cnstr 兜底，这里给足搜索空间
    bound_vel = safety * traj_par.qd_max / sqrt(2);
    bound_acc = safety * traj_par.q2d_max / wk;
    bound_k = min(bound_vel, bound_acc);
    bound_k = max(bound_k, 0.01 * ones(size(bound_k)));
    bound_per_harmonic(:, k) = bound_k;

    idx = (k-1) * 2 * n_dof + (1:2*n_dof);
    ub(idx) = [bound_k; bound_k];
    lb(idx) = -ub(idx);
end

% ab_bound：各关节在所有谐波中的最大系数上界（供参考/打印）
ab_bound = max(bound_per_harmonic, [], 2);

end

%% ===================== 子函数：自动搜索可行初始点 =====================
function [x_best, info] = find_feasible_seed(lb, ub, traj_par, baseQR)
n_opt = numel(lb);
span = ub - lb;% 变量区间长度
center = zeros(n_opt, 1);% 中心点作为初始参考

x_best = center;
[c0, ~] = traj_cnstr(center, traj_par);
info.bestViolation = max(c0);
info.feasibleFound = all(c0 <= 0);
info.bestCost = inf;

if info.feasibleFound
    info.bestCost = traj_cost_lgr(center, traj_par, baseQR);
end

rng(0);
seed_scales = [0.15, 0.30, 0.45, 0.60, 0.85];
trials_per_scale = 12;

for s = seed_scales
    for k = 1:trials_per_scale
        x_try = lb + 0.5 * s * span .* rand(n_opt, 1);
        x_try = min(max(x_try, lb), ub);

        [c_try, ~] = traj_cnstr(x_try, traj_par);
        max_viol = max(c_try);
        if all(c_try <= 0)
            cost_try = traj_cost_lgr(x_try, traj_par, baseQR);
            if ~info.feasibleFound || cost_try < info.bestCost
                x_best = x_try;
                info.bestCost = cost_try;
                info.bestViolation = max_viol;
                info.feasibleFound = true;
            end
        elseif ~info.feasibleFound && max_viol < info.bestViolation
            x_best = x_try;
            info.bestViolation = max_viol;
        end
    end
end

if ~info.feasibleFound
    for alpha = [0.5, 0.25, 0.10, 0.05, 0.01]
        x_try = min(max(alpha * x_best, lb), ub);
        [c_try, ~] = traj_cnstr(x_try, traj_par);
        if all(c_try <= 0)
            x_best = x_try;
            info.bestViolation = max(c_try);
            info.bestCost = traj_cost_lgr(x_try, traj_par, baseQR);
            info.feasibleFound = true;
            break;
        end
    end
end

if ~info.feasibleFound
    warning('find_feasible_seed: 未找到严格可行点，best_violation=%.3e', info.bestViolation);
    x_best = min(max(x_best, lb), ub);
    [c_try, ~] = traj_cnstr(x_best, traj_par);
    info.bestViolation = max(c_try);
    info.bestCost = traj_cost_lgr(x_best, traj_par, baseQR);
end
end
