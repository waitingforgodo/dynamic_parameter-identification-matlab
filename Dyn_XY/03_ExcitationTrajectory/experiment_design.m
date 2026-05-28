% ---------------------------------------------------------------------
% 本脚本用于执行动力学参数辨识的轨迹优化
% ---------------------------------------------------------------------
clear; clc; close all;

% 1. 机器人模型与基础参数加载
addpath(genpath('D:\myWork')); 
robot = importrobot('Marvin.urdf');
% model_type = 'Newton-Euler';
model_type = 'Lagrangian';
[~, baseQR] = GenBaseParameters(robot,model_type);

% 2. 轨迹物理参数设定
traj_par.n_dof = 7;
traj_par.T = 120;                   % 轨迹总时间(s)
traj_par.wf = 2 * pi / traj_par.T;  % 基频率
traj_par.N = 5;                     % 傅里叶级数阶数
traj_par.q0 = deg2rad([0 0 0 -42.5 0 0 0]'); % 初始位置

% 关节限位与运动学约束
traj_par.q_min = deg2rad([-170  -110  -170  -130  -170  -30  -30]');
traj_par.q_max = deg2rad([ 170   110   170    50   170   30   30]');
traj_par.qd_max = [3.0 3.0 3.0 3.0 3.0 3.0 3.0]';
traj_par.q2d_max = [110 100 110 110 110 110 110]';

% 双网格策略：粗网格算目标(提速)，密网格验约束(保安全)
traj_par.t_opt = 0:1.0:traj_par.T;      % 目标函数用粗网格
traj_par.t_cnstr = 0:0.25:traj_par.T;   % 约束检查用密网格
traj_par.t_plot = 0:0.05:traj_par.T;    % 绘图用超密网格
traj_par.t = traj_par.t_cnstr;          % 默认使用密网格
traj_par.t_smp = traj_par.t_cnstr(2) - traj_par.t_cnstr(1);

% 校验初始位置合法性
for i = 1:traj_par.n_dof
    if traj_par.q0(i) < traj_par.q_min(i) || traj_par.q0(i) > traj_par.q_max(i)
        warning('关节 %d 初始位置 q0=%.2f° 超出限位！', i, rad2deg(traj_par.q0(i)));
    end
end

% 3. 优化变量边界估算 (调用子函数)
traj_par.ab_bound = estimate_coeff_bounds(traj_par);
n_opt = 2 * traj_par.n_dof * traj_par.N; % 待优化变量总数

% 按谐波阶数精细化设置上下界 (lb, ub)
lb = zeros(n_opt, 1); ub = zeros(n_opt, 1);
wf = traj_par.wf;
pos_margin = min(traj_par.q_max - traj_par.q0, traj_par.q0 - traj_par.q_min);
pos_margin = max(0.8 * pos_margin, deg2rad(10) * ones(size(pos_margin)));

for k = 1:traj_par.N
    wk = k * wf;
    bound_pos = 0.6 * pos_margin * wk;
    bound_vel = 0.5 * traj_par.qd_max;
    bound_acc = 0.6 * traj_par.q2d_max / wk;
    bound_k = min([bound_pos, bound_vel, bound_acc], [], 2);
    bound_k = max(bound_k, 0.01 * ones(size(bound_k)));
    
    idx = (k-1)*2*traj_par.n_dof + (1:2*traj_par.n_dof);
    ub(idx) = [bound_k; bound_k];
    lb(idx) = -[bound_k; bound_k];
end

% 4. 自动搜索可行初始点 (调用子函数)
fprintf('正在搜索可行的初始种子点...\n');
[x0, seedInfo] = find_feasible_seed(lb, ub, traj_par, baseQR, model_type);
fprintf('寻种完成: 可行=%d, 最佳违反量=%.3e, 最佳代价=%.3f\n\n', ...
    seedInfo.feasibleFound, seedInfo.bestViolation, seedInfo.bestCost);

% 确保并行池已启动
if isempty(gcp('nocreate')), parpool('local'); end


% 5. 执行优化算法
A=[]; b=[]; Aeq=[]; beq=[]; % 无全局线性约束

% 详细配置 GA 参数
optns_ga = optimoptions('ga');
optns_ga.Display = 'iter';                  % 实时播报迭代过程
optns_ga.MaxGenerations = 40;               % 最大进化代数
optns_ga.PopulationSize = 40;               % 种群大小
optns_ga.EliteCount = 8;                    % 精英个体数（保留最好的8个）
optns_ga.FunctionTolerance = 1e-4;          % 目标函数终止容差
optns_ga.InitialPopulationRange = [lb.'; ub.']; % 初始种群生成范围严格限制在 lb/ub 内
optns_ga.PlotFcn = {@gaplotbestf, @gaplotmaxconstr}; % 绘制最佳适应度与最大约束曲线
optns_ga.UseParallel = true;                % 开启并行计算
optns_ga.UseVectorized = false;             % 关闭向量化

% 混合算法精修：GA结束后自动接 patternsearch 进行局部精细挖掘
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
fprintf('开始执行遗传算法(GA)优化...\n');
[x, fval, exitflag, output] = ga(@(x) traj_cost_lgr(x, traj_par, baseQR, model_type), n_opt, ...
    A, b, Aeq, beq, lb, ub, @(x) traj_cnstr(x, traj_par, baseQR, model_type), optns_ga);

% 输出最终结果
[c_final, ~] = traj_cnstr(x, traj_par, baseQR);
fprintf('\n优化结束。Exitflag=%d, 最大约束违反量=%.3e, 最终代价=%.6f\n', ...
    exitflag, max(c_final), fval);

%% ===================== 可选：fmincon 终极局部精调 =====================
refineInfo = struct('used', false, 'accepted', false, 'exitflag', [], 'output', []);

% 检查是否安装了 fmincon，如果有则尝试进行二次精修
if exist('fmincon', 'file') == 2
    try
        fprintf('开始执行 fmincon 终极局部精调...\n');
        optns_fmincon = optimoptions('fmincon', ...
            'Display', 'iter', ...
            'Algorithm', 'sqp', ...               % 序列二次规划算法
            'MaxFunctionEvaluations', 2e4, ...
            'MaxIterations', 120, ...
            'StepTolerance', 1e-4, ...
            'ConstraintTolerance', 1e-6, ...
            'OptimalityTolerance', 1e-4);

        refineInfo.used = true;
        
        % 以 GA 找到的解 x 为起点，进行局部精细挖掘
        [x_ref, fval_ref, exitflag_ref, output_ref] = fmincon( ...
            @(x) traj_cost_lgr(x, traj_par, baseQR, model_type), x, ...
            A, b, Aeq, beq, lb, ub, ...
            @(x) traj_cnstr(x, traj_par, baseQR, model_type), optns_fmincon);

        % 检查精调后的约束与目标函数是否更优
        [c_ref, ~] = traj_cnstr(x_ref, traj_par, baseQR, model_type);
        refineInfo.exitflag = exitflag_ref;
        refineInfo.output = output_ref;

        % 只有当约束满足且目标函数值更小时，才接受精调结果
        if all(c_ref <= 1e-6) && fval_ref < fval
            x = x_ref;
            fval = fval_ref;
            exitflag = exitflag_ref;
            output = output_ref;
            c_final = c_ref;
            refineInfo.accepted = true;
            fprintf('fmincon 精调成功！接受了更优的解。\n');
        else
            fprintf('fmincon 精调未带来明显改善，保留原 GA 结果。\n');
        end
        
    catch ME
        warning('fmincon refinement skipped: %s', ME.message);
    end
end

fprintf('\n最终优化结束。Exitflag=%d, 最大约束违反量=%.3e, 最终代价=%.6f\n', ...
    exitflag, max(c_final), fval);


% 6. 还原轨迹参数并绘图
ab = reshape(x, [2*traj_par.n_dof, traj_par.N]);
a = ab(1:traj_par.n_dof, :); b = ab(traj_par.n_dof+1:2*traj_par.n_dof, :);
c_pol = getPolCoeffs(traj_par.T, a, b, traj_par.wf, traj_par.N, traj_par.q0);
[q, qd, q2d] = mixed_traj(traj_par.t_plot, c_pol, a, b, traj_par.wf, traj_par.N);

jointLabels = arrayfun(@(k) sprintf('q%d', k), 1:traj_par.n_dof, 'UniformOutput', false);
qdLabels = arrayfun(@(k) sprintf('qd%d', k), 1:traj_par.n_dof, 'UniformOutput', false);
q2dLabels = arrayfun(@(k) sprintf('q2d%d', k), 1:traj_par.n_dof, 'UniformOutput', false);

figure('Name', 'Optimized Trajectory')
subplot(3,1,1); plot(traj_par.t_plot, rad2deg(q)); ylabel('$q$', 'interpreter', 'latex'); grid on; legend(jointLabels{:});
subplot(3,1,2); plot(traj_par.t_plot, rad2deg(qd)); ylabel('$\dot{q}$', 'interpreter', 'latex'); grid on; legend(qdLabels{:});
subplot(3,1,3); plot(traj_par.t_plot, rad2deg(q2d)); ylabel('$\ddot{q}$', 'interpreter', 'latex'); grid on; legend(q2dLabels{:}); xlabel('Time (s)');

% 保存结果
currentScriptPath = fileparts(mfilename('fullpath')); 
pathToFolder = fullfile(currentScriptPath, 'results'); 
if ~exist(pathToFolder, 'dir')
    mkdir_status = mkdir(pathToFolder);
    if ~mkdir_status
        error('无法创建文件夹，请检查磁盘权限或路径是否合法: %s', pathToFolder);
    end
end
filename = fullfile(pathToFolder, [model_type, '_7dof_QR.mat']);
save(filename, 'a', 'b', 'c_pol', 'traj_par', 'x', 'fval', 'exitflag', 'output', 'seedInfo', 'refineInfo');
fprintf('轨迹参数已成功保存至: %s\n', filename);


%% ===================== 底部工具子函数区域 =====================

% 基于物理约束估算傅里叶系数边界
function ab_bound = estimate_coeff_bounds(traj_par)
    N = traj_par.N; wf = traj_par.wf;
    % 计算从初始点到各方向限位的最小安全距离
    pos_margin = min(traj_par.q_max - traj_par.q0, traj_par.q0 - traj_par.q_min);
    pos_margin = max(0.8 * pos_margin, deg2rad(10) * ones(size(pos_margin)));

    bound_per_harmonic = zeros(traj_par.n_dof, N);
    for k = 1:N
        wk = k * wf;
        % 根据微分关系推导各阶谐波的最大允许系数
        b_pos = pos_margin * wk / sqrt(N);             % 位置约束
        b_vel = 0.8 * traj_par.qd_max / sqrt(N);       % 速度约束
        b_acc = 0.8 * traj_par.q2d_max / (wk * sqrt(N)); % 加速度约束
        bound_per_harmonic(:, k) = min([b_pos, b_vel, b_acc], [], 2);
    end

    % 取所有谐波中最宽松的那个作为该关节的全局边界
    ab_bound = max(bound_per_harmonic, [], 2);
    ab_bound = max(ab_bound, 0.01 * ones(size(ab_bound))); % 保底最小值
end

% 自动搜索可行初始点
function [x_best, info] = find_feasible_seed(lb, ub, traj_par, baseQR, model_type)
    % 自动搜索一个满足约束条件的可行初始点
    n_opt = numel(lb); span = ub - lb; center = zeros(n_opt, 1);
    x_best = center; 

    % 先测试中心点（零向量）是否可行
    [c0, ~] = traj_cnstr(center, traj_par, baseQR, model_type);
    info.bestViolation = max(c0);
    info.feasibleFound = all(c0 <= 0);
    info.bestCost = inf;
    if info.feasibleFound, info.bestCost = traj_cost_lgr(center, traj_par, baseQR, model_type); end

    rng(0); % 固定随机种子保证可复现
    seed_scales = [0.15, 0.30, 0.45, 0.60]; trials_per_scale = 8;

    % 多尺度随机试探
    for s = seed_scales
        for k = 1:trials_per_scale
            x_try = 0.5 * s * span .* (2 * rand(n_opt, 1) - 1); % 生成随机点
            x_try = min(max(x_try, lb), ub); % 裁剪到边界内
            
            [c_try, ~] = traj_cnstr(x_try, traj_par, baseQR, model_type);
            max_viol = max(c_try);
            
            % 如果找到了完全合法的点，择优录取
            if all(c_try <= 0)
                cost_try = traj_cost_lgr(x_try, traj_par, baseQR, model_type);
                if ~info.feasibleFound || cost_try < info.bestCost
                    x_best = x_try; info.bestCost = cost_try;
                    info.bestViolation = max_viol; info.feasibleFound = true;
                end
            % 如果都没找到合法的，就记录违反程度最小的那个
            elseif ~info.feasibleFound && max_viol < info.bestViolation
                x_best = x_try; info.bestViolation = max_viol;
            end
        end
    end

    % 兜底策略：如果依然没找到可行点，尝试大幅缩小系数幅值
    if ~info.feasibleFound
        for alpha = [0.5, 0.25, 0.10, 0.05, 0.01]
            x_try = alpha * x_best;
            [c_try, ~] = traj_cnstr(x_try, traj_par, baseQR, model_type);
            if all(c_try <= 0)
                x_best = x_try; info.bestViolation = max(c_try);
                info.bestCost = traj_cost_lgr(x_try, traj_par, baseQR, model_type);
                info.feasibleFound = true; break;
            end
        end
    end
end

