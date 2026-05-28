function drvGains = estimate_drive_gains(baseQR, method, path_to_urdf)
% ------------------------------------------------------------------------
% Estimate current-to-torque drive gains using unloaded vs loaded trajectories.
% Supports 7-DOF (requires 7-column currents in CSV and 7-DOF load regressor).
% ------------------------------------------------------------------------
if nargin < 3 || isempty(path_to_urdf)
    path_to_urdf = 'Marvin.urdf';
end

m_load = 2.805;
path_to_unloaded_traj = 'ur-20_02_19_14harm50sec.csv';
path_to_loaded_traj = 'ur-20_02_19_14harm50secLoad.csv';

if isfield(baseQR, 'n_dof') && ~isempty(baseQR.n_dof)
    n_dof = baseQR.n_dof;
else
    n_dof = 7;
end

unloadedTrajectory = parseURData(path_to_unloaded_traj, 195, 4966, n_dof);
unloadedTrajectory = filterData(unloadedTrajectory);

loadedTrajectory = parseURData(path_to_loaded_traj, 308, 5071, n_dof);
loadedTrajectory = filterData(loadedTrajectory);

E1 = baseQR.permutationMatrix(:, 1:baseQR.numberOfBaseParameters);
nFriction = 3 * n_dof;

[Wb_uldd, I_uldd] = buildDriveGainObservation(unloadedTrajectory, baseQR, E1, nFriction);
[Wb_ldd, I_ldd, Wl] = buildDriveGainObservationLoaded(loadedTrajectory, baseQR, E1, nFriction);

Wl_unknown = Wl(:, 1:9);
Wl_known = Wl(:, 10);

if strcmp(method, 'TLS')
    Wb_tls = [I_uldd, -Wb_uldd, zeros(size(I_uldd, 1), size(Wl, 2));
        I_ldd, -Wb_ldd, -Wl_unknown, -Wl_known * m_load];

    [~, ~, V] = svd(Wb_tls, 'econ');
    lmda = 1 / V(end, end);
    pi_tls = lmda * V(:, end);
    drvGains = pi_tls(1:n_dof);
elseif strcmp(method, 'OLS')
    Wb_ls = [I_uldd, -Wb_uldd, zeros(size(I_uldd, 1), size(Wl_unknown, 2));
        I_ldd, -Wb_ldd, -Wl_unknown];

    Yb_ts = [zeros(size(I_uldd, 1), 1); Wl_known * m_load];
    pi_ls = ((Wb_ls' * Wb_ls) \ Wb_ls') * Yb_ts;
    drvGains = pi_ls(1:n_dof);
elseif strcmp(method, 'PC-OLS')
    drvGains = estimateDriveGainsPCOLS(I_uldd, Wb_uldd, I_ldd, Wb_ldd, Wl, ...
        m_load, baseQR, path_to_urdf, n_dof);
else
    error("Chosen method for drive gain estimation does not exist");
end
end


function [Wb, I_diag] = buildDriveGainObservation(traj, baseQR, E1, nFriction)
    n_dof = size(traj.q, 2);
    nSamples = length(traj.t);
    nCols = baseQR.numberOfBaseParameters + nFriction;

    Wb = zeros(nSamples * n_dof, nCols);
    I_diag = zeros(nSamples * n_dof, n_dof);

    for i = 1:nSamples
        q = traj.q(i, :)';
        qd = traj.qd_fltrd(i, :)';
        q2d = traj.q2d_est(i, :)';

        if baseQR.motorDynamicsIncluded
            Yi = regressorWithMotorDynamics(q, qd, q2d);
        else
            Yi = standard_regressor_marvin(q, qd, q2d);
        end
        Ybi = [Yi * E1, frictionRegressor(qd)];

        rows = (i - 1) * n_dof + (1:n_dof);
        Wb(rows, :) = Ybi;
        I_diag(rows, :) = diag(traj.i_fltrd(i, :));
    end
end


function [Wb, I_diag, Wl] = buildDriveGainObservationLoaded(traj, baseQR, E1, nFriction)
    n_dof = size(traj.q, 2);
    nSamples = length(traj.t);
    nCols = baseQR.numberOfBaseParameters + nFriction;

    Wb = zeros(nSamples * n_dof, nCols);
    I_diag = zeros(nSamples * n_dof, n_dof);
    Wl = zeros(nSamples * n_dof, 10);

    for i = 1:nSamples
        q = traj.q(i, :)';
        qd = traj.qd_fltrd(i, :)';
        q2d = traj.q2d_est(i, :)';

        if baseQR.motorDynamicsIncluded
            Yi = regressorWithMotorDynamics(q, qd, q2d);
        else
            Yi = standard_regressor_marvin(q, qd, q2d);
        end
        Ybi = [Yi * E1, frictionRegressor(qd)];
        Yli = load_regressor_UR10E(q, qd, q2d);

        rows = (i - 1) * n_dof + (1:n_dof);
        Wb(rows, :) = Ybi;
        I_diag(rows, :) = diag(traj.i_fltrd(i, :));
        Wl(rows, :) = Yli;
    end
end

% 空载：(K*i_u = W_b*pi)
% 带载：(K*i_l = W_b*pi + W_l*pi_load)
% 相减消去机器人本体参数 π：(K(i_l - i_u) = (W_l*pi_load))
% π_load 里只有质量已知！

function drvGains = estimateDriveGainsPCOLS(I_uldd, Wb_uldd, I_ldd, Wb_ldd, Wl, ...
    m_load, baseQR, path_to_urdf, n_dof)
    bb = baseQR.numberOfBaseParameters;
    nFriction = 3 * n_dof;
    n_link_params = 70;
    if baseQR.motorDynamicsIncluded
        n_std = 11 * n_dof;
    else
        n_std = 10 * n_dof;
    end
    n_dep = n_std - bb;

    drv_gns = sdpvar(n_dof, 1);
    pi_load_unknw = sdpvar(9, 1);
    pi_frctn = sdpvar(nFriction, 1);
    pi_b = sdpvar(bb, 1);
    pi_d = sdpvar(n_dep, 1);

    pii = baseQR.permutationMatrix * ...
        [eye(bb), -baseQR.beta; zeros(n_dep, bb), eye(n_dep)] * [pi_b; pi_d];

    cnstr = [drv_gns(1) > 10];
    mass_indexes = 10:11:n_link_params;
    robot = parse_urdf(path_to_urdf);
    massUpperBound = robot.m(:) * 1.10;

    for i = 1:n_dof
        cnstr = [cnstr, pii(mass_indexes(i)) > 0, ...
            pii(mass_indexes(i)) < massUpperBound(i)]; %#ok<AGROW>
    end

    for i = 1:11:n_link_params
        link_inertia_i = [pii(i), pii(i+1), pii(i+2); ...
            pii(i+1), pii(i+3), pii(i+4); ...
            pii(i+2), pii(i+4), pii(i+5)];
        frst_mmnt_i = vec2skewSymMat(pii(i+6:i+8));
        Di = [link_inertia_i, frst_mmnt_i'; frst_mmnt_i, pii(i+9) * eye(3)];
        cnstr = [cnstr, Di > 0, pii(i+10) > 0]; %#ok<AGROW>
    end

    load_inertia = [pi_load_unknw(1), pi_load_unknw(2), pi_load_unknw(3); ...
        pi_load_unknw(2), pi_load_unknw(4), pi_load_unknw(5); ...
        pi_load_unknw(3), pi_load_unknw(5), pi_load_unknw(6)];
    load_frst_mmnt = vec2skewSymMat(pi_load_unknw(7:9));
    Dl = [load_inertia, load_frst_mmnt'; load_frst_mmnt, m_load * eye(3)];
    cnstr = [cnstr, Dl > 0];

    for i = 1:n_dof
        cnstr = [cnstr, pi_frctn(3*i-2) > 0, pi_frctn(3*i-1) > 0]; %#ok<AGROW>
    end

    t1 = [zeros(size(I_uldd, 1), 1); -Wl(:, end) * m_load];
    t2 = [-I_uldd, Wb_uldd, zeros(size(Wb_uldd, 1), size(Wl, 2) - 1);
        -I_ldd, Wb_ldd, Wl(:, 1:9)];

    obj = norm(t1 - t2 * [drv_gns; pi_b; pi_frctn; pi_load_unknw]);
    optimize(cnstr, obj, sdpsettings('solver', 'sdpt3'));

    drvGains = value(drv_gns);
end
