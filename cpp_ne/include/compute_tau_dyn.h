#pragma once

// Direct torque computation: tau_dyn = Y_base * pi_b (fully expanded).
// Each tau_dyn[i] is a scalar expression of (q, qd, q2d, g).
// Does NOT include friction torque.
void compute_tau_dyn(
    const double q[7],
    const double qd[7],
    const double q2d[7],
    const double g[3],
    double tau_dyn[7]);
