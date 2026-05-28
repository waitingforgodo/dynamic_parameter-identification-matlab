#!/usr/bin/env python3
"""Convert autogen/standard_regressor_marvin.m to C++."""

import re
from pathlib import Path
from typing import List, Optional


def mat_index_to_cpp(line: str) -> str:
    line = re.sub(r"\bin1\((\d+),:\)", lambda m: f"q[{int(m.group(1)) - 1}]", line)
    line = re.sub(r"\bin2\((\d+),:\)", lambda m: f"qd[{int(m.group(1)) - 1}]", line)
    line = re.sub(r"\bin3\((\d+),:\)", lambda m: f"q2d[{int(m.group(1)) - 1}]", line)
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


def mat_expr_to_cpp(expr: str) -> str:
    expr = mat_index_to_cpp(expr)
    expr = replace_elementwise_pow(expr)
    expr = expr.replace(".*", "*")
    expr = expr.replace("./", "/")
    return expr


def convert_line(line: str) -> Optional[str]:
    line = line.strip()
    if not line or line.startswith("%"):
        return None

    if re.match(r"^Y(?:_f)?\s*=\s*reshape", line):
        m = re.search(r"reshape\(\[(.*?)\],7,70\)", line, re.S)
        if not m:
            raise RuntimeError("reshape line not found")
        mt_names = [x.strip() for x in m.group(1).split(",")]
        parts = [
            "    // Column-major reshape to 7x70 (MATLAB compatible)",
            "    int idx = 0;",
            "    const double* blocks[] = {",
        ]
        parts.append("        " + ", ".join(mt_names) + "};")
        parts.append("    const int block_sizes[] = {")
        parts.append("        " + ", ".join(f"static_cast<int>(sizeof({n})/sizeof({n}[0]))" for n in mt_names) + "};")
        parts.extend(
            [
                f"    for (int b = 0; b < {len(mt_names)}; ++b) {{",
                "        for (int k = 0; k < block_sizes[b]; ++k) {",
                "            const int row = idx % 7;",
                "            const int col = idx / 7;",
                "            Y[row * 70 + col] = blocks[b][k];",
                "            ++idx;",
                "        }",
                "    }",
            ]
        )
        return "\n".join(parts)

    m = re.match(r"^(q\d+|q2d\d+|qd\d+|t\d+|et\d+|mt\d+)\s*=\s*(.+);$", line)
    if not m:
        return None

    name, rhs = m.group(1), m.group(2)
    if rhs.startswith("["):
        inner = rhs[1:-1]
        elems = [mat_expr_to_cpp(e.strip()) for e in inner.split(",")]
        return f"    const double {name}[] = {{{', '.join(elems)}}};"
    rhs_cpp = mat_expr_to_cpp(rhs)
    return f"    const double {name} = {rhs_cpp};"


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    src = root / "autogen" / "standard_regressor_marvin.m"
    dst = root / "cpp_le" / "src" / "standard_regressor_marvin.cpp"
    hdr = root / "cpp_le" / "include" / "standard_regressor_marvin.h"

    text = src.read_text(encoding="utf-8")
    reshape_marker = "Y = reshape" if "Y = reshape" in text else "Y_f = reshape"
    start = text.index("q2 = in1(2,:);")
    end = text.rindex(reshape_marker)
    body_lines = text[start : end + text[end:].index(";") + 1].splitlines()

    converted = []  # type: List[str]
    for raw in body_lines:
        out = convert_line(raw)
        if out:
            converted.append(out)

    hdr.write_text(
        """#pragma once

// 7-DOF rigid-body regressor Y(q, qd, qdd), size 7 x 70 (row-major).
void standard_regressor_marvin(
    const double q[7],
    const double dq[7],
    const double ddq[7],
    double Y[7 * 70]);
""",
        encoding="utf-8",
    )

    cpp = (
        "// Auto-converted from autogen/standard_regressor_marvin.m\n"
        '#include "standard_regressor_marvin.h"\n'
        "#include <cmath>\n\n"
        "void standard_regressor_marvin(\n"
        "    const double q[7],\n"
        "    const double dq[7],\n"
        "    const double ddq[7],\n"
        "    double Y[7 * 70])\n"
        "{\n"
        + "\n".join(converted)
        + "\n}\n"
    )
    cpp = re.sub(r"(?<!std::)\bcos\(", "std::cos(", cpp)
    cpp = re.sub(r"(?<!std::)\bsin\(", "std::sin(", cpp)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(cpp, encoding="utf-8")
    print(f"Wrote {dst} ({len(converted)} statements)")


if __name__ == "__main__":
    main()
