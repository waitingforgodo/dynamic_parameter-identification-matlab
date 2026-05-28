function SimuDyn()
clc;
addpath(genpath('D:\myWork'));

thisDir = fileparts(mfilename('fullpath'));
resultDir = fullfile(thisDir, 'result');
modelTypes = {'Newton-Euler', 'Lagrangian'};
runIden = false;   % true=两种方法辨识并保存, false=加载 mat 做双模型验证对比

opts.iden_stride = 5;

if runIden
    runIdentification(resultDir, modelTypes, opts);
else
    runValidation(resultDir, modelTypes);
end

disp('Done');
end

%% ========== 辨识 ==========

function runIdentification(resultDir, modelTypes, opts)
    robot = importrobot('Marvin.urdf');

    if ~exist(resultDir, 'dir')
        mkdir(resultDir);
    end

    for k = 1:numel(modelTypes)
        model_type = modelTypes{k};
        fprintf('\n===== 辨识 [%s] =====\n', model_type);

        [pi_base, baseQR] = GenBaseParameters(robot, model_type);
        [sol, md] = EstimateDynPara(baseQR, model_type, opts);
        sol.pi_base = pi_base;

        matFile = fullfile(resultDir, sprintf('DynIden_%s.mat', model_type));
        save(matFile, 'sol', 'model_type', 'md');
        fprintf('辨识结果已保存: %s\n', matFile);
    end
end

%% ========== 验证对比 ==========

function runValidation(resultDir, modelTypes)
    matNE = fullfile(resultDir, sprintf('DynIden_%s.mat', modelTypes{1}));
    matLG = fullfile(resultDir, sprintf('DynIden_%s.mat', modelTypes{2}));

    if ~isfile(matNE)
        fprintf('错误: 找不到 %s，请设 runIden=true 先运行辨识。\n', matNE);
        return;
    end
    if ~isfile(matLG)
        fprintf('错误: 找不到 %s，请设 runIden=true 先运行辨识。\n', matLG);
        return;
    end

    S_NE = load(matNE);
    S_LG = load(matLG);
    md = S_NE.md;
    fprintf('已加载: %s\n', matNE);
    fprintf('已加载: %s\n', matLG);

    dt = 0.001;
    nSim = size(md.q, 1);
    t = ((1:nSim)' - 1) * dt + 1.002;

    fprintf('\n力矩验证仿真 [%s]...\n', modelTypes{1});
    TauNE = simulateTorque(md, S_NE.sol, modelTypes{1});
    fprintf('\n力矩验证仿真 [%s]...\n', modelTypes{2});
    TauLG = simulateTorque(md, S_LG.sol, modelTypes{2});

    plotDualModelCompare(t, md.tor, TauNE, TauLG, modelTypes);
end

%% ========== 力矩仿真 ==========

function Tau = simulateTorque(md, sol, model_type)
    phi = sol.pi_b;
    pfr = sol.pi_fr;
    E1 = sol.baseQR.permutationMatrix(:, 1:sol.baseQR.numberOfBaseParameters);

    nSim = size(md.q, 1);
    Tau = zeros(nSim, 7);
    updateEvery = max(1, floor(nSim / 200));

    hWait = waitbar(0, sprintf('[%s] 准备计算...', model_type));
    tSim = tic;
    try
        for k = 1:nSim
            if k == 1 || k == nSim || mod(k, updateEvery) == 0
                waitbar(k / nSim, hWait, ...
                    sprintf('[%s] 第 %d / %d 点 (%.1f%%)', model_type, k, nSim, 100 * k / nSim));
            end

            q_rnd = md.q(k, :)';
            dq_rnd = md.dq(k, :)';
            ddq_rnd = md.ddq(k, :)';

            if strcmpi(model_type, 'Newton-Euler')
                Yi = SymForm_NT_Y(q_rnd, dq_rnd, ddq_rnd);
            elseif strcmpi(model_type, 'Lagrangian')
                Yi = SymForm_LG_Y(q_rnd, dq_rnd, ddq_rnd);
            else
                error('Model type error: %s', model_type);
            end

            torD = Yi * E1 * phi;
            torF = SymForm_F(dq_rnd) * pfr;
            Tau(k, :) = (torD + torF)';
        end
    catch ME
        close(hWait);
        rethrow(ME);
    end
    close(hWait);
    fprintf('  %d 点, 耗时 %.1fs\n', nSim, toc(tSim));
end

%% ========== 绘图 ==========

function plotDualModelCompare(t, torMeas, TauNE, TauLG, modelTypes)
    nJoint = size(torMeas, 2);

    for j = 1:nJoint
        figure('Name', sprintf('TorqueCompare_J%d', j));
        plot(t, torMeas(:, j), 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
        hold on;
        plot(t, TauNE(:, j), 'Color', [0 0.45 0.74], 'LineWidth', 1.2);
        plot(t, TauLG(:, j), 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2);

        rmseNE = sqrt(mean((torMeas(:, j) - TauNE(:, j)).^2));
        rmseLG = sqrt(mean((torMeas(:, j) - TauLG(:, j)).^2));
        title(sprintf('Joint %d  RMSE_{NE}=%.4g  RMSE_{LG}=%.4g', j, rmseNE, rmseLG));
        legend('实测转矩', modelTypes{1}, modelTypes{2}, 'Location', 'best');
        xlabel('t (s)');
        ylabel('torque');
        hold off;
    end
end
