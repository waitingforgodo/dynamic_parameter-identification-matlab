# C++ 关节力矩计算（Newton-Euler 回归矩阵）

7-DOF 机械臂力矩前馈计算，基于辨识后的最小参数集。

## 目录结构

```text
cpp_ne/
  include/
    compute_current.h               # 对外接口
      src/

    compute_torque.cpp              # 力矩计算入口 (含摩擦力)
```

## 三级计算接口

| 函数 | 运行时计算量 | 适用场景 |
|------|-------------|---------|
| `compute_joint_torque` | tau_dyn + Y_fr×pi_fr | 参数固定，追求最快速度 |

两个函数都输出含摩擦力的完整力矩。

## 生成 / 更新流程

```bash
# 1. 导出力矩解析式 → C++ (tau_dyn = Y_base * pi_b)
python scripts/convert_tau_dyn_to_cpp.py



```matlab
gen_torque_expression('results/sol.mat')        % → autogen/compute_tau_dyn.m
```

## 使用示例

```cpp
#include "compute_current.h"

double q[7], qd[7], q2d[7], g[3] = {0, 0, 9.8065}, tau[7];
compute_joint_torque(q, qd, q2d, g, PI_FR, tau);
```

## 编译

```bash
cd cpp_ne && mkdir build && cd build
cmake .. && cmake --build .
```
