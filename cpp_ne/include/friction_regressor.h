#pragma once

#include "robot_dyn_config.h"

// Friction regressor: Fv*qd + Fc*sign(qd) + F0 per joint.
// Output Y_fr is 7 x 21 (row-major).
void friction_regressor(const double qd[N_DOF], double Y_fr[N_DOF * N_FRICTION_PARAMS]);
