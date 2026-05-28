#!/usr/bin/env python3
"""Export baseQR_7dof.mat permutation columns (E1) to C++ header."""

from pathlib import Path

import scipy.io as sio


def write_matrix_h(path: Path, name: str, matrix) -> None:
    nr, nc = matrix.shape
    lines = [
        "#pragma once",
        "",
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
    root = Path(__file__).resolve().parents[1]
    mat_path = root / "baseQR_7dof.mat"
    out_dir = root / "cpp_le" / "include"
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

// Exported from baseQR_7dof.mat
static constexpr int N_DOF = {n_dof};
static constexpr int N_STD_PARAMS = {n_std};
static constexpr int N_BASE_PARAMS = {bb};
static constexpr int N_FRICTION_PARAMS = {3 * n_dof};
static constexpr int MOTOR_DYNAMICS_INCLUDED = {motor};
"""
    (out_dir / "robot_dyn_config.h").write_text(config, encoding="utf-8")
    print(f"Exported E1 ({e1.shape[0]}x{e1.shape[1]}) -> {out_dir}")


if __name__ == "__main__":
    main()
