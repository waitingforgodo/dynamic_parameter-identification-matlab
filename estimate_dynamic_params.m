function sol = estimate_dynamic_params(path_to_data, idx, baseQR, method, path_to_urdf, idx_mode)
% -----------------------------------------------------------------------
% 函数功能：从实测轨迹数据中辨识机械臂的基惯性参数与摩擦参数
% 适用对象：7自由度机械臂，支持CSV格式实测数据
% 输入参数：
%   path_to_data   - 实测数据文件路径
%   idx            - 数据截取区间 [起始索引, 结束索引]
%   baseQR         - 最小参数集QR分解结构体（含置换矩阵、基参数数量等）
%   method         - 参数辨识方法:
%                      'OLS'         普通最小二乘
%                      'PC-OLS'      物理一致性约束最小二乘（SDP）
%                      'URDF-REFINE' 以 URDF 原参数为锚点，在保持 OLS 拟合
%                                    精度的前提下做微小修正（推荐用于全参数还原）
%   path_to_urdf   - 机器人URDF模型路径（物理一致性辨识需要）
%   idx_mode       - 数据截取模式 'index' 按点数 / 'time' 按时间
% 输出参数：
%   sol            - 辨识结果结构体，包含基参数、摩擦参数、全参数、标准差等
% -----------------------------------------------------------------------

% 设置URDF文件默认路径
if nargin < 5 || isempty(path_to_urdf)
    path_to_urdf = 'Marvin.urdf';
end

% 设置数据索引模式默认值
if nargin < 6 || isempty(idx_mode)
    idx_mode = 'index';
end

% 第一步：临时读取数据，获取关节数量n_dof
tempData  = parseSimulinkData(path_to_data, idx(1), idx(2), [], idx_mode);
n_dof = tempData.n_dof;

% 第二步：读取指定区间内的辨识轨迹数据
idntfcnTrjctry = parseSimulinkData(path_to_data, idx(1), idx(2), n_dof, idx_mode);

% 第三步：对数据进行滤波（去噪、速度/加速度估计）
idntfcnTrjctry = filterData(idntfcnTrjctry);

% 显示辨识数据的运动范围，用于检查数据是否合理
max_q_deg = rad2deg(max(abs(idntfcnTrjctry.q(:))));
fprintf('辨识数据: max|q|=%.1f deg, max|qd|=%.2f rad/s, max|q2d|=%.1f rad/s^2\n', ...
    max_q_deg, max(abs(idntfcnTrjctry.qd_fltrd(:))), max(abs(idntfcnTrjctry.q2d_est(:))));

% 如果角度超过360度，给出单位异常警告
if max_q_deg > 360
    warning('estimate_dynamic_params: max|q|>360 deg，请检查 parseSimulinkData 单位是否为 deg。');
end

% 第四步：构建观测矩阵Wb和力矩向量Tau
[Tau, Wb] = buildObservationMatrices(idntfcnTrjctry, baseQR);

% ========== 调试输出：观测矩阵状态检查 ==========
fprintf('Wb 矩阵大小: %d x %d\n', size(Wb,1), size(Wb,2));
fprintf('Wb 中是否有 NaN: %d\n', any(isnan(Wb(:))));
fprintf('Wb 秩: %d (应等于列数)\n', rank(Wb));
% ==============================================

% 检测并剔除包含NaN、Inf的异常数据行
bad = any(isnan(Wb), 2) | any(isinf(Wb), 2) | isnan(Tau) | isinf(Tau);
Wb(bad, :) = [];    % 移除异常行对应的回归矩阵行
Tau(bad) = [];      % 移除异常行对应的力矩值

% 初始化输出结果结构体
sol = struct;

% 根据选择的辨识方法执行参数辨识
if strcmp(method, 'OLS')
    % 普通最小二乘辨识（无物理约束）
    [sol.pi_b, sol.pi_fr] = ordinaryLeastSquareEstimation(Tau, Wb, baseQR);
elseif strcmp(method, 'PC-OLS')
    % 带物理约束的最小二乘辨识（质量正定、惯量正定、摩擦正定）
    [sol.pi_b, sol.pi_fr, sol.pi_s] = physicallyConsistentEstimation(...
        Tau, Wb, baseQR, path_to_urdf);
elseif strcmp(method, 'URDF-REFINE')
    % 以 URDF 标称参数为锚点：在 OLS 拟合精度允许范围内，最小化相对 URDF 的偏差
    [sol.pi_b, sol.pi_fr, sol.pi_s, sol.refine_info] = urdfAnchoredEstimation(...
        Tau, Wb, baseQR, path_to_urdf);
else
    % 方法不支持时抛出错误
    error("Chosen method for dynamic parameter estimation does not exist");
end

% ------------------------------------------------------------------
% 辨识结果统计分析：计算参数标准差与相对标准差
% ------------------------------------------------------------------
% 计算辨识误差的无偏方差（残差平方和除以自由度）
% 分子：残差向量的2-范数平方，即实测力矩与预测力矩的误差平方和
% 分母：自由度 = 数据样本数 - 待辨识参数数量，用于无偏估计
sqrd_sgma_e = norm(Tau - Wb * [sol.pi_b; sol.pi_fr], 2)^2 / ...
    (size(Wb, 1) - size(Wb, 2));

% 参数协方差矩阵
% 计算参数的协方差矩阵：方差估计值 × (Wb'*Wb)的伪逆
% 协方差矩阵对角线元素为各参数的方差
Cpi = sqrd_sgma_e * pinv(Wb' * Wb);
% 对协方差矩阵对角线开平方，得到每个辨识参数的【标准差】
sol.std = sqrt(diag(Cpi));    % 各参数的标准差
% 计算【相对标准差】(百分比) = 标准差 / 参数绝对值 × 100
% 用于评估参数辨识的置信度，数值越小表示辨识越可靠
sol.rel_std = 100 * sol.std ./ abs([sol.pi_b; sol.pi_fr]);  % 相对标准差(%)

end


% ======================================================================
% 子函数功能：构建动力学辨识所需的观测矩阵Wb和力矩向量Tau
% 输出 Wb：整体回归矩阵 [基惯性回归项, 摩擦回归项]
% 输出 Tau：实测关节力矩向量
% ======================================================================
function [Tau, Wb] = buildObservationMatrices(idntfcnTrjctry, baseQR)
    % 获取基参数映射矩阵E1：完整参数 → 最小参数集
    E1 = baseQR.permutationMatrix(:, 1:baseQR.numberOfBaseParameters);
    
    n_dof = idntfcnTrjctry.n_dof;           % 关节个数
    nFriction = 3 * n_dof;                  % 摩擦参数个数（每关节3个）
    nSamples = length(idntfcnTrjctry.t);    % 数据采样点总数

    % 初始化观测矩阵Wb：行数=采样点数×关节数，列数=基参数数+摩擦参数数
    Wb = zeros(nSamples * n_dof, baseQR.numberOfBaseParameters + nFriction);
    Tau = zeros(nSamples * n_dof, 1);       % 初始化力矩向量

    % 遍历每一组采样数据，逐点构建回归矩阵
    for i = 1:nSamples
        q = idntfcnTrjctry.q(i, :)';           % 当前时刻关节角度
        qd = idntfcnTrjctry.qd_fltrd(i, :)';   % 滤波后关节速度
        q2d = idntfcnTrjctry.q2d_est(i, :)';   % 估计得到的关节加速度

        % 根据是否包含电机动力学选择对应的回归矩阵函数
        if baseQR.motorDynamicsIncluded
            Yi = regressorWithMotorDynamics(q, qd, q2d);
        else
            Yi = standard_regressor_marvin(q, qd, q2d);
        end
        
        % 计算摩擦回归矩阵
        Yfrctni = frictionRegressor(qd);
        
        % 组合基惯性回归项与摩擦回归项
        Ybi = [Yi * E1, Yfrctni];

        % 计算当前采样点在Wb中的行索引区间
        rows = (i - 1) * n_dof + (1:n_dof);
        Wb(rows, :) = Ybi;              % 填入回归矩阵块
        
        % 填入当前时刻实测滤波后力矩
        Tau(rows) = idntfcnTrjctry.tau_fltrd(i, :)';
    end
end


% ======================================================================
% 子函数功能：普通最小二乘(OLS)参数辨识
% 输出 pib_OLS：基惯性参数
% 输出 pifrctn_OLS：摩擦参数
% ======================================================================
function [pib_OLS, pifrctn_OLS] = ordinaryLeastSquareEstimation(Tau, Wb, baseQR)
    % 最小二乘求解：使用伪逆避免矩阵奇异
    pi_OLS = pinv(Wb) * Tau; 
    
    bb = baseQR.numberOfBaseParameters;    % 基参数数量
    
    pib_OLS     = pi_OLS(1:bb);                  % 提取基惯性参数
    pifrctn_OLS = pi_OLS(bb + 1:end);            % 提取摩擦参数
end


% ======================================================================
% 子函数功能：带物理一致性约束的最小二乘辨识（SDP优化）
% 约束内容：连杆质量>0、惯量矩阵正定、摩擦系数>0
% 输出 pib_SDP：物理可行基参数
% 输出 pifrctn_SDP：物理可行摩擦参数
% 输出 pi_full：还原后的完整惯性参数集
% ======================================================================
function [pib_SDP, pifrctn_SDP, pi_full] = physicallyConsistentEstimation(Tau, Wb, baseQR, path_to_urdf)
    physicalConsistency = 1;    % 使能物理一致性约束

    bb = baseQR.numberOfBaseParameters;        % 基参数数量
    nFriction = size(Wb, 2) - bb;              % 摩擦参数数量
    n_dof = nFriction / 3;                     % 关节数

    % 根据是否包含电机动力学确定标准参数长度
    if baseQR.motorDynamicsIncluded
        n_std = 11 * n_dof;
    else
        n_std = 10 * n_dof;
    end
    n_link_params = n_std;
    n_dep = n_std - bb;    % 依赖参数（冗余参数）的数量

    % 解析URDF模型，获取连杆质量标称值
    robot = parse_urdf(path_to_urdf);
    massValuesURDF = robot.m(1:n_dof);
    disp(massValuesURDF)
    errorRange = 0.010;                        % 质量允许波动范围1%
    massLowerBound = massValuesURDF * (1 - errorRange);
    massUpperBound = massValuesURDF * (1 + errorRange);

    %% 定义SDP优化变量
    pi_frctn = sdpvar(nFriction, 1);        % 摩擦参数变量
    pi_b = sdpvar(bb, 1);                   % 基参数变量
    pi_d = sdpvar(n_dep, 1);                % 依赖参数变量

    % 由基参数+依赖参数还原完整标准惯性参数
    pii = reconstructStandardParams(pi_b, pi_d, baseQR);

    % 构建物理一致性约束（质量>0、惯量正定、摩擦>0）
    cnstr = buildPhysicalConsistencyConstraints(...
        pii, pi_frctn, n_dof, n_link_params, massLowerBound, massUpperBound, physicalConsistency);

    % 优化目标：最小化力矩预测误差
    obj = norm(Tau - Wb * [pi_b; pi_frctn]);

    % 调用SDP求解器进行优化求解
    optimize(cnstr, obj, sdpsettings('solver', 'sdpt3'));

    % 获取优化后的数值解
    pib_SDP = value(pi_b);
    pifrctn_SDP = value(pi_frctn);
    
    % 还原得到完整的连杆惯性参数
    pi_full = reconstructStandardParams(value(pi_b), value(pi_d), baseQR);

end


% ======================================================================
% URDF 锚定辨识：以 URDF 标称参数为中心，在保持动力学精度的前提下
% 输出物理合理、贴近 URDF 的完整惯性参数，用于仿真/可视化/参数还原
% 输入：
%   Tau         - 实测力矩向量
%   Wb          - 观测矩阵（基参数 + 摩擦）
%   baseQR      - 最小参数集结构
%   path_to_urdf- URDF 文件路径
%   opts        - 优化选项（可选）
% 输出：
%   pib         - 优化后的基参数
%   pifrctn     - 优化后的摩擦参数
%   pi_full     - 还原后的完整标准参数
%   info        - 优化信息（收敛性、误差、偏差）
% ======================================================================
function [pib, pifrctn, pi_full, info] = urdfAnchoredEstimation(...
        Tau, Wb, baseQR, path_to_urdf, opts)
    if nargin < 5 || isempty(opts)
        opts = struct();
    end

    fit_rel_tol = getOpt(opts, 'fit_rel_tol', 0.5);
    relBound = getOpt(opts, 'relBound', 0.30);
    mass_rel_bound = getOpt(opts, 'mass_rel_bound', 0.10);
    fit_weight = getOpt(opts, 'fit_weight', 1.0);

    bb = baseQR.numberOfBaseParameters;
    nFriction = size(Wb, 2) - bb;
    n_dof = nFriction / 3;
    n_link_params = baseQR.n_std_params;
    n_dep = n_link_params - bb;

    pi_std_urdf = getUrdfStandardParams(path_to_urdf, n_dof, baseQR);
    [pi_b_urdf, pi_d_urdf] = splitBaseDepFromStandard(pi_std_urdf, baseQR);
    robot = parse_urdf(path_to_urdf);
    massLowerBound = robot.m(1:n_dof)' * (1 - mass_rel_bound);
    massUpperBound = robot.m(1:n_dof)' * (1 + mass_rel_bound);

    pi_ols = pinv(Wb) * Tau;
    fit_ols = norm(Tau - Wb * pi_ols);
    fit_limit_sq = ((1 + fit_rel_tol) * fit_ols)^2;

    % 正规方程：||Wb*pi-Tau||^2 = pi'*H*pi - 2*g'*pi + Tau'*Tau，仅 70x70，避免 68 万行 SOC
    H = Wb' * Wb;
    g = Wb' * Tau;
    tauNormSq = Tau' * Tau;
    fit_rhs = fit_limit_sq - tauNormSq;

    delta_b = sdpvar(bb, 1);
    delta_d = sdpvar(n_dep, 1);
    pi_fr = sdpvar(nFriction, 1);
    pi_all = [pi_b_urdf + delta_b; pi_fr];

    pi_b = pi_b_urdf + delta_b;
    pi_d = pi_d_urdf + delta_d;
    pii = reconstructStandardParams(pi_b, pi_d, baseQR);

    cnstr = buildPhysicalConsistencyConstraints(...
        pii, pi_fr, n_dof, n_link_params, massLowerBound, massUpperBound, true);
    % cnstr = [cnstr, buildUrdfDeviationBounds(pii, pi_std_urdf, relBound)]; %#ok<AGROW>
    cnstr = [cnstr, pi_all' * H * pi_all - 2 * g' * pi_all <= fit_rhs]; %#ok<AGROW>

    pi_fr_ols = pi_ols(bb + 1:end);
    fr_upper = max(abs(pi_fr_ols) * 5, 1.0);
    cnstr = [cnstr, pi_fr <= fr_upper]; %#ok<AGROW>

    param_scale = max(abs(pi_std_urdf), 1e-6);
    anchor_obj = norm((pii - pi_std_urdf) ./ param_scale);
    obj = anchor_obj + fit_weight * (pi_all' * H * pi_all - 2 * g' * pi_all) / max(tauNormSq, 1e-12);

    sol = solveUrdfRefineSdp(cnstr, obj);

    if sol.problem ~= 0
        fprintf('urdfAnchoredEstimation: SDP 失败 (code=%d, %s)，尝试 fmincon 备选...\n', ...
            sol.problem, yalmiperror(sol.problem));
        [pib, pifrctn, pi_d_opt, ok] = urdfRefineFminconFallback(...
            Tau, Wb, baseQR, pi_std_urdf, pi_b_urdf, pi_d_urdf, ...
            fit_ols, fit_rel_tol, relBound, mass_rel_bound, massLowerBound, massUpperBound);
        if ~ok
            fprintf('urdfAnchoredEstimation: fmincon 亦失败，返回 OLS 解。\n');
            pib = pi_ols(1:bb);
            pifrctn = pi_ols(bb + 1:end);
            pi_full = pi_std_urdf;
            info = struct('converged', false, 'problem', sol.problem, ...
                'fit_ols', fit_ols / sqrt(numel(Tau)), ...
                'fit_refined', norm(Tau - Wb * pi_ols) / sqrt(numel(Tau)));
            return;
        end
        pi_full = reconstructStandardParams(pib, pi_d_opt, baseQR);
        fit_refined = norm(Tau - Wb * [pib; pifrctn]);
        rel_dev = abs(pi_full - pi_std_urdf) ./ param_scale;
        info = packRefineInfo(true, fit_ols, fit_refined, fit_rel_tol, relBound, ...
            rel_dev, pi_std_urdf, pi_d_opt, 'fmincon', numel(Tau));
        printRefineSummary(info, fit_rel_tol);
        return;
    end

    pib = value(pi_b);
    pifrctn = value(pi_fr);
    pi_d_opt = value(pi_d);
    pi_full = reconstructStandardParams(pib, pi_d_opt, baseQR);

    fit_refined = norm(Tau - Wb * [pib; pifrctn]);
    rel_dev = abs(pi_full - pi_std_urdf) ./ param_scale;
    info = packRefineInfo(true, fit_ols, fit_refined, fit_rel_tol, relBound, ...
        rel_dev, pi_std_urdf, pi_d_opt, 'sdp', numel(Tau));
    printRefineSummary(info, fit_rel_tol);
end


function sol = solveUrdfRefineSdp(cnstr, obj)
    solvers = {'sdpt3', 'sedumi'};
    sol = struct('problem', -3);
    for k = 1:numel(solvers)
        sol = optimize(cnstr, obj, sdpsettings('solver', solvers{k}, 'verbose', 0));
        if sol.problem == 0
            return;
        end
    end
end


function [pib, pifrctn, pi_d, ok] = urdfRefineFminconFallback(...
        Tau, Wb, baseQR, pi_std_urdf, pi_b_urdf, pi_d_urdf, ...
        fit_ols, fit_rel_tol, relBound, mass_rel_bound, massLowerBound, massUpperBound)
    bb = baseQR.numberOfBaseParameters;
    nFriction = size(Wb, 2) - bb;
    n_dof = nFriction / 3;
    n_link_params = baseQR.n_std_params;
    n_dep = n_link_params - bb;
    fit_limit = (1 + fit_rel_tol) * fit_ols;
    param_scale = max(abs(pi_std_urdf), 1e-6);

    pi_ols = pinv(Wb) * Tau;
    x0 = [zeros(bb + n_dep, 1); max(pi_ols(bb + 1:end), 0.01)];
    lb = [-relBound * abs([pi_b_urdf; pi_d_urdf]); zeros(nFriction, 1)];
    ub = [relBound * max(abs([pi_b_urdf; pi_d_urdf]), 1e-3); inf(nFriction, 1)];

    function f = objFun(x)
        pi_b = pi_b_urdf + x(1:bb);
        pi_d = pi_d_urdf + x(bb + 1:bb + n_dep);
        pi_fr = x(bb + n_dep + 1:end);
        pii = reconstructStandardParams(pi_b, pi_d, baseQR);
        fit_term = norm(Wb * [pi_b; pi_fr] - Tau);
        anchor_term = norm((pii - pi_std_urdf) ./ param_scale);
        f = anchor_term + fit_term / max(norm(Tau), 1e-12);
    end

    function [c, ceq] = nlcon(x)
        pi_b = pi_b_urdf + x(1:bb);
        pi_d = pi_d_urdf + x(bb + 1:bb + n_dep);
        pi_fr = x(bb + n_dep + 1:end);
        pii = reconstructStandardParams(pi_b, pi_d, baseQR);
        ceq = [];
        c = [norm(Wb * [pi_b; pi_fr] - Tau) - fit_limit];
        mass_idx = 10:10:n_link_params;
        mass_idx = mass_idx(1:n_dof);
        for i = 1:n_dof
            c(end + 1) = massLowerBound(i) - pii(mass_idx(i)); %#ok<AGROW>
            c(end + 1) = pii(mass_idx(i)) - massUpperBound(i); %#ok<AGROW>
        end
        for i = 1:n_dof
            c(end + 1) = -pi_fr(3*i-2); %#ok<AGROW>
            c(end + 1) = -pi_fr(3*i-1); %#ok<AGROW>
        end
        for i = 1:numel(pi_std_urdf)
            delta_max = relBound * param_scale(i);
            c(end + 1) = abs(pii(i) - pi_std_urdf(i)) - delta_max; %#ok<AGROW>
        end
    end

    opts = optimoptions('fmincon', 'Display', 'off', 'MaxIterations', 500, ...
        'Algorithm', 'sqp', 'ConstraintTolerance', 1e-6);
    [x, ~, exitflag] = fmincon(@objFun, x0, [], [], [], [], lb, ub, @nlcon, opts);
    ok = exitflag > 0;
    pib = pi_b_urdf + x(1:bb);
    pi_d = pi_d_urdf + x(bb + 1:bb + n_dep);
    pifrctn = x(bb + n_dep + 1:end);
end


function info = packRefineInfo(converged, fit_ols, fit_refined, fit_rel_tol, relBound, ...
        rel_dev, pi_std_urdf, pi_d, solverName, nTau)
    info = struct();
    info.converged = converged;
    info.solver = solverName;
    info.fit_ols = fit_ols / sqrt(nTau);
    info.fit_refined = fit_refined / sqrt(nTau);
    info.fit_rel_tol = fit_rel_tol;
    info.relBound = relBound;
    info.max_rel_dev_std = max(rel_dev);
    info.mean_rel_dev_std = mean(rel_dev);
    info.pi_std_urdf = pi_std_urdf;
    info.pi_d = pi_d;
end


function printRefineSummary(info, fit_rel_tol)
    fit_limit = (1 + fit_rel_tol) * info.fit_ols;
    fprintf('URDF-REFINE [%s]: OLS RMSE=%.4f, refined RMSE=%.4f (limit %.4f)\n', ...
        info.solver, info.fit_ols, info.fit_refined, fit_limit);
    fprintf('URDF-REFINE: 全参数相对 URDF 偏差 max=%.2f%%, mean=%.2f%%\n', ...
        100 * info.max_rel_dev_std, 100 * info.mean_rel_dev_std);
end

% ======================================================================
% 从URDF获取标准参数向量
% ======================================================================
function pi_std = getUrdfStandardParams(path_to_urdf, n_dof, baseQR)
    robot = parse_urdf(path_to_urdf);
    % 如果包含电机动力学，补充电机参数
    if baseQR.motorDynamicsIncluded
        robot.pi = [robot.pi; zeros(1, n_dof)];
    end
    pi_std = robot.pi(:);
end


% ======================================================================
% 将完整标准参数分解为 基参数 + 依赖参数
% ======================================================================
function [pi_b, pi_d] = splitBaseDepFromStandard(pi_std, baseQR)
    bb = baseQR.numberOfBaseParameters;
    E = baseQR.permutationMatrix;
    beta = baseQR.beta;
    perm = E' * pi_std;
    pi_d = perm(bb + 1:end);
    pi_b = perm(1:bb) + beta * pi_d;
end


% ======================================================================
% 函数功能：由【基参数 pi_b】 + 【依赖参数 pi_d】
% 还原得到机器人完整的标准惯性参数集 pi_std（10个/连杆）
% 这是最小参数集 → 完整参数集的核心映射公式
% 输入：
%   pi_b   - 辨识得到的基参数（最小参数集）
%   pi_d   - 依赖参数（冗余参数）
%   baseQR - 最小参数分解结构体
% 输出：
%   pi_std - 完整标准惯性参数向量（可直接用于URDF/多体仿真）
% ======================================================================
function pi_std = reconstructStandardParams(pi_b, pi_d, baseQR)
    % 基参数（可辨识参数）的数量
    bb = baseQR.numberOfBaseParameters;
    
    % 依赖参数（冗余参数）的数量
    n_dep = numel(pi_d);
    
    % 核心公式：最小参数 → 完整惯性参数
    % 变换步骤：
    % pi_std = 置换矩阵 * 变换矩阵 * [基参数; 依赖参数]
    pi_std = baseQR.permutationMatrix * ...
        [eye(bb), -baseQR.beta; zeros(n_dep, bb), eye(n_dep)] * [pi_b; pi_d];
end


% ======================================================================
% 构建物理一致性约束：质量>0、惯量正定、摩擦>0
% ======================================================================
function cnstr = buildPhysicalConsistencyConstraints(...
        pii, pi_frctn, n_dof, n_link_params, massLowerBound, massUpperBound, physicalConsistency)
    cnstr = [];
    % 质量参数索引：每10个参数第10个为质量
    mass_indexes = 10:10:n_link_params;
    mass_indexes = mass_indexes(1:n_dof);

    % 约束质量在 URDF 标称值允许波动范围内
    for i = 1:n_dof
        cnstr = [cnstr, pii(mass_indexes(i)) > massLowerBound(i), ...
            pii(mass_indexes(i)) < massUpperBound(i)]; %#ok<AGROW>
    end

    % 惯性张量物理约束（Dyad矩阵正定）
    if physicalConsistency
        for i = 1:10:n_link_params
            link_inertia_i = [pii(i), pii(i+1), pii(i+2); ...
                pii(i+1), pii(i+3), pii(i+4); ...
                pii(i+2), pii(i+4), pii(i+5)];
            frst_mmnt_i = pii(i+6:i+8);
            Di = [0.5 * trace(link_inertia_i) * eye(3) - link_inertia_i, ...
                frst_mmnt_i; frst_mmnt_i', pii(i+9)];
            cnstr = [cnstr, Di > 1e-6 * eye(4)]; %#ok<AGROW>
        end
    else
        % 简化版惯性约束
        for i = 1:10:n_link_params
            link_inertia_i = [pii(i), pii(i+1), pii(i+2); ...
                pii(i+1), pii(i+3), pii(i+4); ...
                pii(i+2), pii(i+4), pii(i+5)];
            frst_mmnt_i = vec2skewSymMat(pii(i+6:i+8));
            Di = [link_inertia_i, frst_mmnt_i'; frst_mmnt_i, pii(i+9) * eye(3)];
            cnstr = [cnstr, Di >= 0]; %#ok<AGROW>
        end
    end

    % 摩擦参数约束：粘性、库仑摩擦 > 0
    for i = 1:n_dof
        cnstr = [cnstr, pi_frctn(3*i-2) > 0, pi_frctn(3*i-1) > 0]; %#ok<AGROW>
    end
end


% ======================================================================
% 构建相对于URDF的参数偏差约束
% ======================================================================
function cnstr = buildUrdfDeviationBounds(pii, pi_std_urdf, relBound)
    cnstr = [];
    scale = max(abs(pi_std_urdf), 1e-6);
    for i = 1:numel(pi_std_urdf)
        delta_max = relBound * scale(i);
        cnstr = [cnstr, ...
            pii(i) >= pi_std_urdf(i) - delta_max, ...
            pii(i) <= pi_std_urdf(i) + delta_max]; %#ok<AGROW>
    end
end


% ======================================================================
% 获取可选参数，无则返回默认值
% ======================================================================
function val = getOpt(opts, name, defaultVal)
    if isfield(opts, name) && ~isempty(opts.(name))
        val = opts.(name);
    else
        val = defaultVal;
    end
end