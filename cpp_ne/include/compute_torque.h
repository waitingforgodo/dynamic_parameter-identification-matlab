#pragma once

#include "robot_dyn_config.h"

// tau = [Y*E1, Y_fr] * [pi_b; pi_fr]
// Same method as validate_dynamic_params.m
void compute_joint_torque(
    const double q[N_DOF],
    const double qd[N_DOF],
    const double q2d[N_DOF],
    const double g[N_GRAVITY_INPUTS],
    const double pi_b[N_BASE_PARAMS],
    const double pi_fr[N_FRICTION_PARAMS],
    double tau[N_DOF]);
