function [pi_lgr_base, baseQR] = base_params_qr(includeMotorDynamics, cacheFile)
% ----------------------------------------------------------------------
% QR decomposition to find base (minimum) parameters from regressor.
% Optimized: pre-allocation, fewer samples, optional caching.
%
% Usage:
%   [~, baseQR] = base_params_qr(0);                   % compute fresh
%   [~, baseQR] = base_params_qr(0, 'baseQR.mat');    % load if exists
% ----------------------------------------------------------------------
% Seed the random number generator based on the current time
% rng('shuffle');
rng(0);
n_dof = 7;

% --- Check cache ---
if nargin >= 2 && ~isempty(cacheFile) && isfile(cacheFile)
    tmp = load(cacheFile, 'baseQR', 'pi_lgr_base');
    if tmp.baseQR.motorDynamicsIncluded == includeMotorDynamics
        baseQR = tmp.baseQR;
        pi_lgr_base = tmp.pi_lgr_base;
        fprintf('Loaded baseQR from cache: %s\n', cacheFile);
        return;
    end
end

% rng(0);  % fixed seed for reproducibility
robot=importrobot('Marvin.urdf');
q_min = [robot.Bodies{1,1}.Joint.PositionLimits(1) robot.Bodies{1,2}.Joint.PositionLimits(1) robot.Bodies{1,3}.Joint.PositionLimits(1) ...
         robot.Bodies{1,4}.Joint.PositionLimits(1) robot.Bodies{1,5}.Joint.PositionLimits(1) robot.Bodies{1,6}.Joint.PositionLimits(1) ...
         robot.Bodies{1,7}.Joint.PositionLimits(1)]';
q_max = [robot.Bodies{1,1}.Joint.PositionLimits(2) robot.Bodies{1,2}.Joint.PositionLimits(2) robot.Bodies{1,3}.Joint.PositionLimits(2) ...
         robot.Bodies{1,4}.Joint.PositionLimits(2) robot.Bodies{1,5}.Joint.PositionLimits(2) robot.Bodies{1,6}.Joint.PositionLimits(2) ...
         robot.Bodies{1,7}.Joint.PositionLimits(2)]';
qd_max  = 3 *q_max;
q2d_max = 6 * q_max;
% --- Joint limits ---
% q_min = -pi * ones(n_dof, 1);
% q_max =  pi * ones(n_dof, 1);
% qd_max  = 3 * pi * ones(n_dof, 1);
% q2d_max = 6 * pi * ones(n_dof, 1);

% --- Symbolic parameters (only needed for pi_lgr_base output) ---
m   = sym('m%d',   [n_dof,1], 'real');
hx  = sym('h%d_x', [n_dof,1], 'real');
hy  = sym('h%d_y', [n_dof,1], 'real');
hz  = sym('h%d_z', [n_dof,1], 'real');
ixx = sym('i%d_xx',[n_dof,1], 'real');
ixy = sym('i%d_xy',[n_dof,1], 'real');
ixz = sym('i%d_xz',[n_dof,1], 'real');
iyy = sym('i%d_yy',[n_dof,1], 'real');
iyz = sym('i%d_yz',[n_dof,1], 'real');
izz = sym('i%d_zz',[n_dof,1], 'real');
im  = sym('im%d',  [n_dof,1], 'real');

if includeMotorDynamics
    nLnkPrms = 11;
    pi_lgr_sym = [ixx, ixy, ixz, iyy, iyz, izz, hx, hy, hz, m, im]';
else
    nLnkPrms = 10;
    pi_lgr_sym = [ixx, ixy, ixz, iyy, iyz, izz, hx, hy, hz, m]';
    % pi_lgr_sym = [m,hx, hy, hz,ixx, ixy, ixz, iyy, iyz, izz]';
end
pi_lgr_sym = pi_lgr_sym(:);
nLnks = n_dof;

% --- Build observation matrix W (pre-allocated) ---
nSamples = 25;
nCols = nLnkPrms * n_dof;
W = zeros(nSamples * n_dof, nCols);

fprintf('Building observation matrix (%d samples)...\n', nSamples);
for i = 1:nSamples
    q_rnd   = q_min + (q_max - q_min) .* rand(n_dof, 1);
    qd_rnd  = -qd_max + 2 * qd_max .* rand(n_dof, 1);
    q2d_rnd = -q2d_max + 2 * q2d_max .* rand(n_dof, 1);

    if includeMotorDynamics
        Y = regressorWithMotorDynamics(q_rnd, qd_rnd, q2d_rnd);
    else
        g=[0 0 9.8065]';
        Y = standard_regressor_marvin(q_rnd, qd_rnd, q2d_rnd,g);
    end
    rows = (i - 1) * n_dof + (1:n_dof);
    W(rows, :) = Y;
end

% --- QR decomposition with column pivoting ---
[~, R, E] = qr(W);

% Determine rank from R diagonal
tol = max(size(W)) * eps(norm(R(1,1)));
bb = sum(abs(diag(R)) > 1e-8);
% bb=rank(W);
R1 = R(1:bb, 1:bb);
R2 = R(1:bb, bb+1:end);
beta = R1 \ R2;
beta(abs(beta) < sqrt(eps)) = 0;

% Verify relation W2 = W1 * beta
W1 = W * E(:, 1:bb);
W2 = W * E(:, bb+1:end);
relErr = norm(W2 - W1 * beta, 'fro') / (norm(W1, 'fro') + eps);
assert(relErr < 1e-6, 'QR relation check failed (relErr=%.2e)', relErr);

% --- Symbolic base parameters ---
pi1 = E(:, 1:bb)' * pi_lgr_sym;
pi2 = E(:, bb+1:end)' * pi_lgr_sym;
pi_lgr_base = pi1 + beta * pi2;
% disp(pi_lgr_base);

% --- Output structure ---
baseQR = struct;
baseQR.numberOfBaseParameters = bb;
baseQR.permutationMatrix = E;
baseQR.beta = beta;
baseQR.motorDynamicsIncluded = includeMotorDynamics;
baseQR.n_dof = nLnks;
baseQR.n_std_params = nLnkPrms * nLnks;

fprintf('Found %d base parameters (rank=%d)\n', bb, bb);

% --- Save cache if requested ---
if nargin >= 2 && ~isempty(cacheFile)
    save(cacheFile, 'baseQR', 'pi_lgr_base');
    fprintf('Saved baseQR to cache: %s\n', cacheFile);
end

end
