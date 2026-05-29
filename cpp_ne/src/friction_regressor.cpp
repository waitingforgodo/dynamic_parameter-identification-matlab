#include "kdl/friction_regressor.h"

#include <cmath>

void friction_regressor(const double qd[N_DOF], double Y_fr[N_DOF * N_FRICTION_PARAMS]) {
    for (int i = 0; i < N_DOF * N_FRICTION_PARAMS; ++i) {
        Y_fr[i] = 0.0;
    }
    for (int i = 0; i < N_DOF; ++i) {
        const double sgn = (qd[i] > 0.0) ? 1.0 : ((qd[i] < 0.0) ? -1.0 : 0.0);
        const int base_col = 3 * i;
        Y_fr[i * N_FRICTION_PARAMS + base_col + 0] = qd[i];
        Y_fr[i * N_FRICTION_PARAMS + base_col + 1] = sgn;
        Y_fr[i * N_FRICTION_PARAMS + base_col + 2] = 1.0;
    }
}
