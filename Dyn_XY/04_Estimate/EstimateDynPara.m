function [sol, md] = EstimateDynPara(baseQR, model_type, opts)
% 最小二乘动力学参数辨识
% opts.iden_stride     采样步长，默认 10（每 10 点取 1 点，约 10 倍加速）
% opts.iden_max_samples 最大采样点数，默认 inf

if nargin < 3, opts = struct(); end
cfg.iden_stride = 10;
cfg.iden_max_samples = inf;
fn = fieldnames(opts);
for k = 1:numel(fn), cfg.(fn{k}) = opts.(fn{k}); end

[md] = MotionDataProcess();
[nRow, ~] = size(md.q);

idxSample = 1:cfg.iden_stride:nRow;
if isfinite(cfg.iden_max_samples) && numel(idxSample) > cfg.iden_max_samples
    idxSample = idxSample(1:cfg.iden_max_samples);
end

nBase = baseQR.numberOfBaseParameters;
nParam = nBase + 21;
E1 = baseQR.permutationMatrix(:, 1:nBase);

A = zeros(nParam, nParam);
b = zeros(nParam, 1);
tauSq = 0;
nObs = 0;

tBuild = tic;
for k = 1:numel(idxSample)
    i = idxSample(k);
    q_rnd = md.q(i, :)';
    dq_rnd = md.dq(i, :)';
    ddq_rnd = md.ddq(i, :)';

    if strcmpi(model_type, 'Newton-Euler')
        Yi = SymForm_NT_Y(q_rnd, dq_rnd, ddq_rnd);
    elseif strcmpi(model_type, 'Lagrangian')
        Yi = SymForm_LG_Y(q_rnd, dq_rnd, ddq_rnd);
    else
        error('Model type error: %s', model_type);
    end

    Ybi = [Yi * E1, SymForm_F(dq_rnd)];
    torList = md.tor(i, :)';

    A = A + Ybi' * Ybi;
    b = b + Ybi' * torList;
    tauSq = tauSq + torList' * torList;
    nObs = nObs + 7;
end
tBuild = toc(tBuild);

pi_OLS = A \ b;

sol = struct();
sol.pi_b = pi_OLS(1:nBase);
sol.pi_fr = pi_OLS(nBase+1:end);
sol.baseQR = baseQR;
sol.iden_cfg = cfg;
sol.iden_idxSample = idxSample;

ssRes = max(tauSq - pi_OLS' * b, 0);
rmse = sqrt(ssRes / nObs);
relErr = sqrt(ssRes) / (sqrt(tauSq) + eps);

fprintf('辨识采样: stride=%d, %d/%d 点, 构建 %.1fs\n', ...
    cfg.iden_stride, numel(idxSample), nRow, tBuild);
fprintf('辨识完成 [%s]: obs=%d, pi_b=%d, pi_fr=%d, RMSE=%.4g, relErr=%.2e\n', ...
    model_type, nObs, nBase, numel(sol.pi_fr), rmse, relErr);

end
