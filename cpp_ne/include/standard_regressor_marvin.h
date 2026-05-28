#pragma once

// SymForm_NT_YG: Y(q, qd, qdd, g), size 7 x 70 (row-major).
// g is base-frame gravity linear acceleration [m/s^2] (e.g. {0, 0, 9.8065}).
void standard_regressor_marvin(
    const double q[7],
    const double qd[7],
    const double q2d[7],
    const double g[3],
    double Y[7 * 70]);
