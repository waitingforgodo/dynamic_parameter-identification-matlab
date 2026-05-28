#!/usr/bin/env python3
"""Convert autogen/standard_regressor_marvin.m (SymForm_NT_YG) to C++.

Supports the 4-argument regressor Y(q, qd, qdd, g) from GenRegNewtonEulerGravity.
Keeps scripts/convert_regressor_to_cpp.py unchanged (legacy 3-arg format).
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List, Optional, Tuple


def mat_index_to_cpp(line: str) -> str:
    line = re.sub(r"\bin1\((\d+),:\)", lambda m: f"q[{int(m.group(1)) - 1}]", line)
    line = re.sub(r"\bin2\((\d+),:\)", lambda m: f"qd[{int(m.group(1)) - 1}]", line)
    line = re.sub(r"\bin3\((\d+),:\)", lambda m: f"q2d[{int(m.group(1)) - 1}]", line)
    line = re.sub(r"\bin4\((\d+),:\)", lambda m: f"g[{int(m.group(1)) - 1}]", line)
    return line


def replace_elementwise_pow(expr: str) -> str:
    while True:
        m = re.search(r"([\w]+)\.\^(\d+(?:\.\d+)?)", expr)
        if not m:
            break
        base, exp = m.group(1), m.group(2)
        if exp in {"2", "2.0"}:
            repl = f"({base})*({base})"
        elif exp in {"3", "3.0"}:
            repl = f"({base})*({base})*({base})"
        else:
            repl = f"std::pow({base}, {exp})"
        expr = expr[: m.start()] + repl + expr[m.end() :]
    return expr


def unpack_symvars_to_cpp(expr: str) -> str:
    """Map MATLAB unpacked names (q1, dq1, ddq1, g1) to C array indices."""
    for i in range(7, 0, -1):
        expr = re.sub(rf"\bddq{i}\b", f"q2d[{i - 1}]", expr)
    for i in range(7, 0, -1):
        expr = re.sub(rf"\bdq{i}\b", f"qd[{i - 1}]", expr)
    for i in range(7, 0, -1):
        expr = re.sub(rf"\bq{i}\b", f"q[{i - 1}]", expr)
    for i in range(3, 0, -1):
        expr = re.sub(rf"\bg{i}\b", f"g[{i - 1}]", expr)
    return expr


def mat_expr_to_cpp(expr: str) -> str:
    expr = mat_index_to_cpp(expr)
    expr = unpack_symvars_to_cpp(expr)
    expr = replace_elementwise_pow(expr)
    expr = expr.replace(".*", "*")
    expr = expr.replace("./", "/")
    return expr


def parse_reshape_dims(text: str) -> Tuple[int, int, List[str]]:
    m = re.search(
        r"Y\s*=\s*reshape\(\[(.*?)\],\s*(\d+)\s*,\s*(\d+)\s*\)",
        text,
        re.S,
    )
    if not m:
        raise RuntimeError("Y = reshape([...], n_rows, n_cols) not found")
    mt_names = [x.strip() for x in m.group(1).split(",") if x.strip()]
    return int(m.group(2)), int(m.group(3)), mt_names


def convert_reshape_line(n_rows: int, n_cols: int, mt_names: List[str]) -> str:
    parts = [
        f"    // Column-major reshape to {n_rows}x{n_cols} (MATLAB compatible)",
        "    int idx = 0;",
        "    const double* blocks[] = {",
        "        " + ", ".join(mt_names) + "};",
        "    const int block_sizes[] = {",
        "        "
        + ", ".join(f"static_cast<int>(sizeof({n})/sizeof({n}[0]))" for n in mt_names)
        + "};",
        f"    for (int b = 0; b < {len(mt_names)}; ++b) {{",
        "        for (int k = 0; k < block_sizes[b]; ++k) {",
        f"            const int row = idx % {n_rows};",
        f"            const int col = idx / {n_rows};",
        f"            Y[row * {n_cols} + col] = blocks[b][k];",
        "            ++idx;",
        "        }",
        "    }",
    ]
    return "\n".join(parts)


def convert_line(line: str, n_rows: int, n_cols: int) -> Optional[str]:
    line = line.strip()
    if not line or line.startswith("%"):
        return None

    if re.match(r"^Y\s*=\s*reshape", line):
        _, _, mt_names = parse_reshape_dims(line if "reshape" in line else f"Y = {line}")
        return convert_reshape_line(n_rows, n_cols, mt_names)

    m = re.match(
        r"^(q\d+|dq\d+|ddq\d+|g\d+|t\d+|et\d+|mt\d+)\s*=\s*(.+);$",
        line,
    )
    if not m:
        return None

    name, rhs = m.group(1), m.group(2)
    if rhs.startswith("["):
        inner = rhs[1:-1]
        elems = [mat_expr_to_cpp(e.strip()) for e in inner.split(",")]
        return f"    const double {name}[] = {{{', '.join(elems)}}};"
    rhs_cpp = mat_expr_to_cpp(rhs)
    return f"    const double {name} = {rhs_cpp};"


def find_body_range(text: str) -> Tuple[int, int, int, int]:
    n_rows, n_cols, _ = parse_reshape_dims(text)
    start_markers = ("t2 = cos(q1);", "t2=cos(q1);")
    start = -1
    for marker in start_markers:
        if marker in text:
            start = text.index(marker)
            break
    if start < 0:
        raise RuntimeError("Could not locate start of regressor body (expected t2 = cos(q1))")
    end = text.rindex("Y = reshape")
    end = text.index(";", end) + 1
    return start, end, n_rows, n_cols


def write_header(path: Path, n_rows: int, n_cols: int) -> None:
    path.write_text(
        f"""#pragma once

// SymForm_NT_YG: Y(q, qd, qdd, g), size {n_rows} x {n_cols} (row-major).
// g is base-frame gravity linear acceleration [m/s^2] (e.g. {{0, 0, 9.8065}}).
void standard_regressor_marvin(
    const double q[7],
    const double qd[7],
    const double q2d[7],
    const double g[3],
    double Y[{n_rows} * {n_cols}]);
""",
        encoding="utf-8",
    )


def write_cpp(
    path: Path,
    converted: List[str],
    n_rows: int,
    n_cols: int,
    src_name: str,
) -> None:
    body = "\n".join(converted)
    cpp = (
        f"// Auto-converted from {src_name} (SymForm_NT_YG, gravity input g[3])\n"
        '#include "standard_regressor_marvin.h"\n'
        "#include <cmath>\n\n"
        "void standard_regressor_marvin(\n"
        "    const double q[7],\n"
        "    const double qd[7],\n"
        "    const double q2d[7],\n"
        "    const double g[3],\n"
        f"    double Y[{n_rows} * {n_cols}])\n"
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
        help="MATLAB regressor file (default: autogen/standard_regressor_marvin.m)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="C++ tree with include/ and src/ (default: cpp_ne)",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    src = args.src or (root / "autogen" / "standard_regressor_marvin.m")
    out_dir = args.out_dir or (root / "cpp_ne")
    dst = out_dir / "src" / "standard_regressor_marvin.cpp"
    hdr = out_dir / "include" / "standard_regressor_marvin.h"

    text = src.read_text(encoding="utf-8")
    if "in4" not in text and "SymForm_NT_YG" not in text:
        raise RuntimeError(
            f"{src} does not look like SymForm_NT_YG (missing in4). "
            "Use convert_regressor_to_cpp.py for the legacy 3-arg regressor."
        )

    start, end, n_rows, n_cols = find_body_range(text)
    body_lines = text[start:end].splitlines()

    converted: List[str] = []
    for raw in body_lines:
        out = convert_line(raw, n_rows, n_cols)
        if out:
            converted.append(out)

    write_header(hdr, n_rows, n_cols)
    write_cpp(dst, converted, n_rows, n_cols, src.relative_to(root).as_posix())
    print(f"Wrote {hdr}")
    print(f"Wrote {dst} ({len(converted)} statements, Y is {n_rows}x{n_cols})")


if __name__ == "__main__":
    main()
