function out = getPolCoeffs(T, a, b, wf, N, q0)
% -----------------------------------------------------------------------
% Compute 5th-order polynomial coefficients for endpoint continuity.
% Engineering improvement: cache the constant linear system matrix.
% -----------------------------------------------------------------------

n = size(q0, 1);
persistent Ac_cache T_cache n_cache

if isempty(Ac_cache) || isempty(T_cache) || isempty(n_cache) || T_cache ~= T || n_cache ~= n
    I = eye(n);
    O = zeros(n);
    Ac_cache = [I,  O,  O,  O,  O,  O; ...
                O,  I,  O,  O,  O,  O; ...
                O,  O,  2*I, O,  O,  O; ...
                I,  T*I, T^2*I, T^3*I, T^4*I, T^5*I; ...
                O,  I,  2*T*I, 3*T^2*I, 4*T^3*I, 5*T^4*I; ...
                O,  O,  2*I, 6*T*I, 12*T^2*I, 20*T^3*I];
    T_cache = T;
    n_cache = n;
end

qh0 = -sum(b ./ ((1:N) * wf), 2);
qdh0 = sum(a, 2);
q2dh0 = sum(b .* ((1:N) * wf), 2);

[qhT, qdhT, q2dhT] = fourier_series_traj(T, zeros(n, 1), a, b, wf, N);
rhs = [q0 - qh0; -qdh0; -q2dh0; q0 - qhT; -qdhT; -q2dhT];

c = Ac_cache \ rhs;
out = reshape(c, [n, 6]);
end
