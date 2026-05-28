function out = traj_cost_lgr(opt_vars, traj_par, baseQR)
% 轨迹优化代价（越小越好）：-logdet + 利用率/运动/对称惩罚

opt_vars = opt_vars(:);
N = traj_par.N;
wf = traj_par.wf;
T = traj_par.T;
n = size(traj_par.q0, 1);

if isfield(traj_par, 't_opt') && ~isempty(traj_par.t_opt)
    t_logdet = reshape(traj_par.t_opt, 1, []);
else
    t_logdet = reshape(traj_par.t, 1, []);
end
if isfield(traj_par, 't_cnstr') && ~isempty(traj_par.t_cnstr)
    t_kin = reshape(traj_par.t_cnstr, 1, []);
else
    t_kin = t_logdet;
end

ab = reshape(opt_vars, [2*n, N]);
a = ab(1:n, :);
b = ab(n+1:2*n, :);
c_pol = getPolCoeffs(T, a, b, wf, N, traj_par.q0);

% logdet 用 t_opt；峰值利用率用更密的 t_cnstr，避免漏掉速度/加速度峰值
[q_k, qd_k, q2d_k] = mixed_traj(t_kin, c_pol, a, b, wf, N);
[q, qd, q2d] = mixed_traj(t_logdet, c_pol, a, b, wf, N);

viol = normalizedLimitViolation(q_k, qd_k, q2d_k, traj_par);
viol_pos = max(viol, 0);
motion_level = mean(std(qd_k, 0, 2) ./ max(traj_par.qd_max(:), 1e-6));

if any(viol_pos > 0)
    out = 1e6 * (1 + sum(viol_pos.^2)) + 1e3 / (motion_level + 1e-6);
    return;
end
% 
% E1 = baseQR.permutationMatrix(:, 1:baseQR.numberOfBaseParameters);
% nFriction = 3 * n;
% W = zeros(numel(t_logdet) * n, baseQR.numberOfBaseParameters + nFriction);
% 
% for i = 1:numel(t_logdet)
%     if baseQR.motorDynamicsIncluded
%         Yi = regressorWithMotorDynamics(q(:, i), qd(:, i), q2d(:, i));
%     else
%         Yi = standard_regressor_marvin(q(:, i), qd(:, i), q2d(:, i));
%     end
%     W((i-1)*n + (1:n), :) = [Yi * E1, frictionRegressor(qd(:, i))];
% end
% 
% if any(~isfinite(W(:)))
%     out = 1e9;
%     return;
% end
% 
% col_norms = sqrt(sum(W.^2, 1));
% col_norms(col_norms < 1e-10) = 1;
% Wn = W ./ col_norms;
% F = Wn' * Wn + 1e-9 * eye(size(Wn, 2));
% 
% if any(~isfinite(F(:)))
%     out = 1e9;
%     return;
% end
% 
% [R, p] = chol(F);
% if p ~= 0
%     out = 1e8 + trace(F);
%     return;
% end
% logdetF = 2 * sum(log(diag(R) + eps));

% --- 惩罚项 ---
w_motion = 15.0;
if isfield(traj_par, 'w_motion'), w_motion = traj_par.w_motion; end
motion_penalty = w_motion / (motion_level + 1e-3);

w_sym = 6.0;
if isfield(traj_par, 'w_symmetry'), w_sym = traj_par.w_symmetry; end
q_above = max(q_k - traj_par.q0, [], 2);
q_below = max(traj_par.q0 - q_k, [], 2);
symmetry_ratio = min(q_above, q_below) ./ (max(q_above, q_below) + 1e-6);
symmetry_penalty = w_sym * sum(1 - symmetry_ratio);

target_qd  = 0.65;
target_q2d = 0.60;
w_u_qd  = 80.0;
w_u_q2d = 60.0;
if isfield(traj_par, 'target_qd_util'),  target_qd  = traj_par.target_qd_util;  end
if isfield(traj_par, 'target_q2d_util'), target_q2d = traj_par.target_q2d_util; end
if isfield(traj_par, 'w_util_qd'),  w_u_qd  = traj_par.w_util_qd;  end
if isfield(traj_par, 'w_util_q2d'), w_u_q2d = traj_par.w_util_q2d; end

qd_peak  = max(abs(qd_k), [], 2) ./ max(traj_par.qd_max(:), 1e-6);
q2d_peak = max(abs(q2d_k), [], 2) ./ max(traj_par.q2d_max(:), 1e-6);
% 平方惩罚：峰值距目标越远，代价增长越快
util_penalty = w_u_qd  * sum(max(0, target_qd  - qd_peak).^2) + ...
               w_u_q2d * sum(max(0, target_q2d - q2d_peak).^2);

regularization = 1e-6 * (opt_vars' * opt_vars) / numel(opt_vars);
% out = -logdetF + motion_penalty + symmetry_penalty + util_penalty + regularization;
out =motion_penalty + symmetry_penalty + util_penalty + regularization;
end

function viol = normalizedLimitViolation(q, qd, q2d, traj_par)
q_span  = max(traj_par.q_max(:) - traj_par.q_min(:), 1e-6);
qd_lim  = max(traj_par.qd_max(:), 1e-6);
q2d_lim = max(traj_par.q2d_max(:), 1e-6);

viol = [
    (traj_par.q_min(:) - min(q, [], 2)) ./ q_span;
    (max(q, [], 2) - traj_par.q_max(:)) ./ q_span;
    (max(abs(qd), [], 2) - traj_par.qd_max(:)) ./ qd_lim;
    (max(abs(q2d), [], 2) - traj_par.q2d_max(:)) ./ q2d_lim];
end
