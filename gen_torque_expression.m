function gen_torque_expression(solFile, cacheFile)
%GEN_TORQUE_EXPRESSION 将基回归矩阵与辨识参数相乘，生成力矩符号表达式
%
%   将 Y_base(q, dq, ddq, g) * pi_b 在符号层面直接展开，
%   生成 tau_dyn(q, dq, ddq, g) 的 MATLAB 函数文件。
%
%   这样 C++ 中惯性力矩计算变为直接的解析表达式，
%   无需任何矩阵乘法（既不需要 E1，也不需要 pi_b 的点乘循环）。
%
%   用法:
%     gen_torque_expression()
%     gen_torque_expression('results/sol.mat', 'baseQR_7dof.mat')
%
%   输出: autogen/compute_tau_dyn.m
%         tau = compute_tau_dyn(q, dq, ddq, g)
%         返回 7x1 向量，为惯性/重力/科氏力矩（不含摩擦力）

if nargin < 1 || isempty(solFile)
    solFile = fullfile('results', 'sol.mat');
end
if nargin < 2 || isempty(cacheFile)
    cacheFile = 'baseQR_7dof.mat';
end

fprintf('=== 生成力矩解析表达式 tau_dyn = Y_base * pi_b ===\n');

%% 1. 加载辨识结果 sol
fprintf('加载辨识结果: %s\n', solFile);
tmp = load(solFile, 'sol');
sol = tmp.sol;
pi_b = sol.pi_b(:);
bb = numel(pi_b);
fprintf('  基参数数量: %d\n', bb);

%% 2. 加载 baseQR 获取 E1 矩阵
fprintf('加载 baseQR: %s\n', cacheFile);
tmp = load(cacheFile, 'baseQR');
baseQR = tmp.baseQR;
E1 = baseQR.permutationMatrix(:, 1:baseQR.numberOfBaseParameters);
n_dof = baseQR.n_dof;
n_std = baseQR.n_std_params;

assert(bb == baseQR.numberOfBaseParameters, ...
    'sol.pi_b 长度(%d) 与 baseQR.numberOfBaseParameters(%d) 不匹配', ...
    bb, baseQR.numberOfBaseParameters);

%% 3. 构建符号变量
q_sym = sym('q%d', [n_dof, 1], 'real');
dq_sym = sym('dq%d', [n_dof, 1], 'real');
ddq_sym = sym('ddq%d', [n_dof, 1], 'real');
g_sym = sym('g', [3, 1], 'real');

%% 4. 计算符号标准回归矩阵
fprintf('计算符号标准回归矩阵 Y_std (%dx%d)...\n', n_dof, n_std);
Y_std = standard_regressor_marvin(q_sym, dq_sym, ddq_sym, g_sym);

%% 5. 符号乘法 Y_base = Y_std * E1
fprintf('计算 Y_base = Y_std * E1 (%dx%d)...\n', n_dof, bb);
Y_base = Y_std * E1;

%% 6. 将数值 pi_b 代入，得到 tau_dyn = Y_base * pi_b (7x1 符号向量)
fprintf('计算 tau_dyn = Y_base * pi_b (%dx1)...\n', n_dof);
tau_dyn = Y_base * pi_b;

%% 7. 化简
% fprintf('化简表达式...\n');
% tau_dyn = expand(tau_dyn);

%% 8. 生成 MATLAB 函数文件
repoRoot = fileparts(mfilename('fullpath'));
autogenDir = fullfile(repoRoot, 'autogen');
if ~exist(autogenDir, 'dir')
    mkdir(autogenDir);
end
outFile = fullfile(autogenDir, 'compute_tau_dyn');
fprintf('生成函数文件: %s.m\n', outFile);
matlabFunction(tau_dyn, 'File', outFile, 'Vars', {q_sym, dq_sym, ddq_sym, g_sym});

%% 9. 同时保存 pi_b 信息供后续验证
infoFile = fullfile(autogenDir, 'compute_tau_dyn_info.mat');
save(infoFile, 'pi_b', 'solFile', 'cacheFile');

fprintf('=== 完成! ===\n');
fprintf('  调用: tau = compute_tau_dyn(q, dq, ddq, g)\n');
fprintf('  tau 大小: %d x 1 (惯性+重力+科氏力矩，不含摩擦)\n', n_dof);
fprintf('  C++ 中每个关节力矩为一个直接的解析表达式，无需矩阵运算\n');

end
