#!/usr/bin/env python3
"""Export baseQR_7dof.mat (SymForm_NT_YG / NE+gravity QR) to C++ headers.

Reads the root-level baseQR_7dof.mat produced after regenerating with g-aware
regressor. Keeps scripts/export_baseqr_cpp.py unchanged (legacy cache paths).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import scipy.io as sio


def write_matrix_h(path: Path, name: str, matrix) -> None:
    nr, nc = matrix.shape
    lines = [
        "#pragma once",
        "",
        f"// Exported from baseQR_7dof.mat (NE regressor with gravity input g)",
        f"static constexpr int {name}_ROWS = {nr};",
        f"static constexpr int {name}_COLS = {nc};",
        f"static constexpr double {name}[{name}_ROWS][{name}_COLS] = {{",
    ]
    for r in range(nr):
        row = ", ".join(f"{matrix[r, c]:.17e}" for c in range(nc))
        lines.append(f"    {{{row}}},")
    lines.append("};")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mat",
        type=Path,
        default=None,
        help="baseQR .mat file (default: baseQR_7dof.mat in repo root)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Directory for base_qr_E1.h and robot_dyn_config.h (default: cpp_ne/include)",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    mat_path = args.mat or (root / "baseQR_7dof.mat")
    out_dir = args.out_dir or (root / "cpp_ne" / "include")
    out_dir.mkdir(parents=True, exist_ok=True)

    data = sio.loadmat(str(mat_path))
    base_qr = data["baseQR"][0, 0]
    e_full = base_qr["permutationMatrix"]
    bb = int(base_qr["numberOfBaseParameters"][0, 0])
    n_dof = int(base_qr["n_dof"][0, 0])
    n_std = int(base_qr["n_std_params"][0, 0])
    motor = int(base_qr["motorDynamicsIncluded"][0, 0])

    e1 = e_full[:, :bb]
    write_matrix_h(out_dir / "base_qr_E1.h", "E1", e1)

    config = f"""#pragma once

// Exported from {mat_path.name} (SymForm_NT_YG + base_params_qr)
// Gravity g is a regressor INPUT (3-vector), not a column in Y ({n_std} std params).
static constexpr int N_DOF = {n_dof};
static constexpr int N_STD_PARAMS = {n_std};
static constexpr int N_BASE_PARAMS = {bb};
static constexpr int N_FRICTION_PARAMS = {3 * n_dof};
static constexpr int MOTOR_DYNAMICS_INCLUDED = {motor};
static constexpr int N_GRAVITY_INPUTS = 3;
"""
    (out_dir / "robot_dyn_config.h").write_text(config, encoding="utf-8")
    print(f"Exported E1 ({e1.shape[0]}x{e1.shape[1]}) -> {out_dir}")
    print(f"N_BASE_PARAMS={bb}, N_STD_PARAMS={n_std}")


if __name__ == "__main__":
    main()
