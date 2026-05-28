# C++ 关节力矩计算（Newton-Euler 回归矩阵）

参考 `validate_dynamic_params.m` 中的方法，在 C++ 中实现：

```text
tau = [Y(q,qd,qdd) * E1,  Y_fr(qd)] * [pi_b; pi_fr]
```

其中 `E1 = baseQR.permutationMatrix(:, 1:numberOfBaseParameters)`，来自 `baseQR_7dof.mat`。

## 目录结构

```text
cpp_ne/
  include/
    robot_dyn_config.h      # 维度常量（7 DOF, 49 基参数, 21 摩擦参数）
    base_qr_E1.h            # 从 baseQR_7dof.mat 导出的 E1 (70×49)
    standard_regressor_marvin.h
    friction_regressor.h
    compute_torque.h        # 对外接口：compute_joint_torque(...)
    dynamic_params.h        # pi_b, pi_fr（辨识后导出）
  src/
    standard_regressor_marvin.cpp   # 由 autogen/standard_regressor_marvin.m 自动转换
    friction_regressor.cpp
    compute_torque.cpp
```

## 生成 / 更新数据

```bash
# 1. 从 baseQR_7dof.mat 导出 E1
python scripts/export_baseqr_cpp.py

# 2. 从 autogen/standard_regressor_marvin.m 生成 C++ 回归矩阵
python scripts/convert_regressor_to_cpp.py

# 3. 在 MATLAB 中辨识并导出 pi_b、pi_fr
export_dynamic_params_cpp(sol);
```

## 在工程中使用

```cpp
#include "compute_torque.h"
#include "dynamic_params.h"

double q[7], qd[7], q2d[7], tau[7];
// ... 填入关节位置、速度、加速度（弧度制）...
compute_joint_torque(q, qd, q2d, PI_B, PI_FR, tau);
```

## 编译

```bash
cd cpp_ne && mkdir build && cd build
cmake ..
cmake --build .
```

生成静态库 `robot_dyn_ne`，可在你的工程中链接使用。

## 核心 API

`compute_joint_torque(q, qd, q2d, pi_b, pi_fr, tau)` — 输入 7 关节 `q/qd/q2d` 与辨识参数，输出 7 维预测力矩 (Nm)。
