#pragma once

// Exported from baseQR_7dof.mat (SymForm_NT_YG + base_params_qr)
// Gravity g is a regressor INPUT (3-vector), not a column in Y (70 std params).
static constexpr int N_DOF = 7;
static constexpr int N_STD_PARAMS = 70;
static constexpr int N_BASE_PARAMS = 49;
static constexpr int N_FRICTION_PARAMS = 21;
static constexpr int MOTOR_DYNAMICS_INCLUDED = 0;
static constexpr int N_GRAVITY_INPUTS = 3;
