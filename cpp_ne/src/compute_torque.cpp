#include "compute_torque.h"

#include "base_qr_E1.h"
#include "friction_regressor.h"
#include "standard_regressor_marvin.h"

void compute_joint_torque(
    const double q[N_DOF],
    const double qd[N_DOF],
    const double q2d[N_DOF],
    const double g[N_GRAVITY_INPUTS],
    const double pi_b[N_BASE_PARAMS],
    const double pi_fr[N_FRICTION_PARAMS],
    double tau[N_DOF]) {
    double Y_std[N_DOF * N_STD_PARAMS];
    double Y_fr[N_DOF * N_FRICTION_PARAMS];
    standard_regressor_marvin(q, qd, q2d, g, Y_std);
    friction_regressor(qd, Y_fr);

    for (int i = 0; i < N_DOF; ++i) {
        double sum = 0.0;
        for (int j = 0; j < N_BASE_PARAMS; ++j) {
            double yb = 0.0;
            for (int k = 0; k < N_STD_PARAMS; ++k) {
                yb += Y_std[i * N_STD_PARAMS + k] * E1[k][j];
            }
            sum += yb * pi_b[j];
        }
        for (int j = 0; j < N_FRICTION_PARAMS; ++j) {
            sum += Y_fr[i * N_FRICTION_PARAMS + j] * pi_fr[j];
        }
        tau[i] = sum;
    }
}
