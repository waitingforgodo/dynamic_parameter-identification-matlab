#!/usr/bin/env python3
"""Convert autogen/compute_tau_dyn.m to C++.

This is the fully-expanded torque expression: tau = Y_base * pi_b.
Each joint torque is a direct scalar expression of (q, qd, q2d, g),
with no matrix operations at runtime.

Handles matlabFunction output format where:
  mt1 = [expr1; expr2; expr3; expr4];
  mt2 = [expr5; expr6; expr7];
  tau_dyn = [mt1; mt2];
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def unpack_symvars_to_cpp(expr: str) -> str:
    for i in range(7, 0, -1):
        expr = re.sub(rf"\bddq{i}\b", f"q2d[{i - 1}]", expr)
    for i in range(7, 0, -1):
        expr = re.sub(rf"\bdq{i}\b", f"qd[{i - 1}]", expr)
    for i in range(7, 0, -1):
        expr = re.sub(rf"\bq{i}\b", f"q[{i - 1}]", expr)
    for i in range(3, 0, -1):
        expr = re.sub(rf"\bg{i}\b", f"g[{i - 1}]", expr)
    return expr


def mat_index_to_cpp(line: str) -> str:
    line = re.sub(r"\bin1\((\d+),:\)", lambda m: f"q[{int(m.group(1)) - 1}]", line)
    line = re.sub(r"\bin2\((\d+),:\)", lambda m: f"qd[{int(m.group(1)) - 1}]", line)
    line = re.sub(r"\bin3\((\d+),:\)", lambda m: f"q2d[{int(m.group(1)) - 1}]", line)
    line = re.sub(r"\bin4\((\d+),:\)", lambda m: f"g[{int(m.group(1)) - 1}]", line)
    return line


def replace_elementwise_pow(expr: str) -> str:
    while True:
        m = re.search(r"([\w\[\]]+)\.\^(\d+(?:\.\d+)?)", expr)
        if not m:
            break
        base, exp_str = m.group(1), m.group(2)
        if exp_str in {"2", "2.0"}:
            repl = f"({base})*({base})"
        elif exp_str in {"3", "3.0"}:
            repl = f"({base})*({base})*({base})"
        else:
            repl = f"std::pow({base}, {exp_str})"
        expr = expr[: m.start()] + repl + expr[m.end() :]
    return expr


def mat_expr_to_cpp(expr: str) -> str:
    expr = mat_index_to_cpp(expr)
    expr = unpack_symvars_to_cpp(expr)
    expr = replace_elementwise_pow(expr)
    expr = expr.replace(".*", "*")
    expr = expr.replace("./", "/")
    return expr


def parse_mt_elements(rhs: str) -> List[str]:
    """Parse [expr1;expr2;...] into list of expressions."""
    inner = rhs.strip()
    if inner.startswith("[") and inner.endswith("]"):
        inner = inner[1:-1]
    return [e.strip() for e in inner.split(";") if e.strip()]


def convert_file(text: str) -> Tuple[List[str], int]:
    """Convert the MATLAB function body to C++ statements."""
    converted: List[str] = []
    mt_elements: Dict[str, List[str]] = {}

    lines = text.splitlines()

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("%") or stripped.startswith("function"):
            continue
        if stripped == "end":
            continue

        # Skip variable unpacking lines
        if re.match(r"^(q\d+|dq\d+|ddq\d+|g\d+)\s*=\s*in\d+\(", stripped):
            continue

        # mt array assignment: mt1 = [expr1;expr2;...]
        m = re.match(r"^(mt\d+)\s*=\s*(\[.+\]);$", stripped)
        if m:
            name = m.group(1)
            rhs = m.group(2)
            elements = parse_mt_elements(rhs)
            mt_elements[name] = elements
            for i, elem in enumerate(elements):
                elem_cpp = mat_expr_to_cpp(elem)
                converted.append(f"    const double {name}_{i} = {elem_cpp};")
            continue

        # tau_dyn = [mt1;mt2] or tau_dyn = [mt1,mt2]
        m = re.match(r"^tau_dyn\s*=\s*\[(.*?)\];$", stripped)
        if m:
            inner = m.group(1)
            parts = [p.strip() for p in re.split(r"[;,]", inner) if p.strip()]
            idx = 0
            for part in parts:
                if part in mt_elements:
                    for i in range(len(mt_elements[part])):
                        converted.append(
                            f"    tau_dyn[{idx}] = {part}_{i};"
                        )
                        idx += 1
                else:
                    elem_cpp = mat_expr_to_cpp(part)
                    converted.append(f"    tau_dyn[{idx}] = {elem_cpp};")
                    idx += 1
            return converted, idx

        # Temp variable assignment (t2, et1, etc.)
        m = re.match(r"^(t\d+|et\d+)\s*=\s*(.+);$", stripped)
        if m:
            name, rhs = m.group(1), m.group(2)
            rhs_cpp = mat_expr_to_cpp(rhs)
            converted.append(f"    const double {name} = {rhs_cpp};")
            continue

    return converted, 0


def write_header(path: Path, n_dof: int) -> None:
    path.write_text(
        f"""#pragma once

// Direct torque computation: tau_dyn = Y_base * pi_b (fully expanded).
// Each tau_dyn[i] is a scalar expression of (q, qd, q2d, g).
// Does NOT include friction torque.
void compute_tau_dyn(
    const double q[{n_dof}],
    const double qd[{n_dof}],
    const double q2d[{n_dof}],
    const double g[3],
    double tau_dyn[{n_dof}]);
""",
        encoding="utf-8",
    )


def write_cpp(
    path: Path,
    converted: List[str],
    n_dof: int,
    src_name: str,
) -> None:
    body = "\n".join(converted)
    cpp = (
        f"// Auto-converted from {src_name}\n"
        f"// tau_dyn = Y_base * pi_b (fully expanded, no matrix ops at runtime)\n"
        '#include "compute_tau_dyn.h"\n'
        "#include <cmath>\n\n"
        "void compute_tau_dyn(\n"
        f"    const double q[{n_dof}],\n"
        f"    const double qd[{n_dof}],\n"
        f"    const double q2d[{n_dof}],\n"
        "    const double g[3],\n"
        f"    double tau_dyn[{n_dof}])\n"
        "{\n"
        + body
        + "\n}\n"
    )
    cpp = re.sub(r"(?<!std::)\bcos\(", "std::cos(", cpp)
    cpp = re.sub(r"(?<!std::)\bsin\(", "std::sin(", cpp)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(cpp, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--src",
        type=Path,
        default=None,
        help="MATLAB tau_dyn file (default: autogen/compute_tau_dyn.m)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="C++ tree with include/ and src/ (default: cpp_ne)",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    src = args.src or (root / "autogen" / "compute_tau_dyn.m")
    out_dir = args.out_dir or (root / "cpp_ne")
    dst = out_dir / "src" / "compute_tau_dyn.cpp"
    hdr = out_dir / "include" / "compute_tau_dyn.h"

    if not src.exists():
        raise FileNotFoundError(
            f"{src} not found. Run gen_torque_expression.m in MATLAB first."
        )

    text = src.read_text(encoding="utf-8")
    converted, n_dof = convert_file(text)

    if not converted:
        raise RuntimeError("No statements converted. Check input file format.")
    if n_dof == 0:
        raise RuntimeError("Could not determine n_dof from tau_dyn assignment.")

    write_header(hdr, n_dof)
    write_cpp(dst, converted, n_dof, src.relative_to(root).as_posix())
    print(f"Wrote {hdr}")
    print(f"Wrote {dst} ({len(converted)} statements, n_dof={n_dof})")


if __name__ == "__main__":
    main()
