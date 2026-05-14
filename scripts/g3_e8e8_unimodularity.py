from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

from golden_octonions import PHI, q_basis_icosian


Z_BASIS_LABELS = ["1", "i", "h", "g", "phi", "phi*i", "phi*h", "phi*g"]


def z_trace_polar_matrix():
    r_basis = q_basis_icosian()
    z_basis = r_basis + [b.scale(PHI) for b in r_basis]
    matrix = []
    for x in z_basis:
        row = []
        for y in z_basis:
            value = (x.inner(y) * 2).trace_q()
            if value.denominator != 1:
                raise ValueError(f"non-integral trace-polar entry: {value}")
            row.append(int(value))
        matrix.append(row)
    return matrix


def block_diag(a, b):
    n = len(a)
    m = len(b)
    out = [[0 for _ in range(n + m)] for _ in range(n + m)]
    for i in range(n):
        for j in range(n):
            out[i][j] = a[i][j]
    for i in range(m):
        for j in range(m):
            out[n + i][n + j] = b[i][j]
    return out


def leading_principal_minors(matrix):
    return [int(sp.Matrix([row[:i] for row in matrix[:i]]).det()) for i in range(1, len(matrix) + 1)]


def smith_diagonal(matrix):
    snf = smith_normal_form(sp.Matrix(matrix))
    return [int(abs(snf[i, i])) for i in range(min(snf.shape)) if snf[i, i] != 0]


def standard_e8_gram():
    return [
        [2, -1, 0, 0, 0, 0, 0, 0],
        [-1, 2, -1, 0, 0, 0, 0, 0],
        [0, -1, 2, -1, 0, 0, 0, 0],
        [0, 0, -1, 2, -1, 0, 0, 0],
        [0, 0, 0, -1, 2, -1, 0, -1],
        [0, 0, 0, 0, -1, 2, -1, 0],
        [0, 0, 0, 0, 0, -1, 2, 0],
        [0, 0, 0, 0, -1, 0, 0, 2],
    ]


def build_certificate():
    icosian = z_trace_polar_matrix()
    g0 = block_diag(icosian, icosian)
    e8 = standard_e8_gram()
    icosian_det = int(sp.Matrix(icosian).det())
    g0_det = int(sp.Matrix(g0).det())
    e8_det = int(sp.Matrix(e8).det())
    icosian_even = all(icosian[i][i] % 2 == 0 for i in range(8))
    g0_even = all(g0[i][i] % 2 == 0 for i in range(16))
    icosian_positive = all(x > 0 for x in leading_principal_minors(icosian))
    g0_positive = all(x > 0 for x in leading_principal_minors(g0))
    return {
        "certificate": "e8e8_unimodularity_certificate",
        "basis_convention": {
            "R_basis": ["1", "i", "h=(1+i+j+k)/2", "g=(-1+(phi-1)i-phi*j)/2"],
            "Z_basis": Z_BASIS_LABELS,
            "trace_polar_form": "Tr_{K/Q}(2*<x,y>_K)",
        },
        "icosian_gram_matrix": icosian,
        "icosian_gram_determinant": icosian_det,
        "icosian_gram_smith_diagonal": smith_diagonal(icosian),
        "icosian_form_even": icosian_even,
        "icosian_form_positive_definite_by_sylvester": icosian_positive,
        "g0_gram_matrix": g0,
        "g0_gram_determinant": g0_det,
        "g0_gram_smith_diagonal": smith_diagonal(g0),
        "g0_form_even": g0_even,
        "g0_form_positive_definite_by_sylvester": g0_positive,
        "standard_e8_gram_matrix": e8,
        "standard_e8_determinant": e8_det,
        "explicit_isometry_to_E8": {
            "exists": False,
            "matrix": None,
            "reason": "Determinant mismatch: the trace-polar icosian Z-lattice has determinant 625, while E8 has determinant 1.",
        },
        "decision_gate": {
            "expected_det_1_confirmed": icosian_det == 1 and g0_det == 1,
            "normalization_risk_triggered": icosian_det != 1 or g0_det != 1,
            "modified_discriminant_group_analysis_required": icosian_det > 1 and int(icosian_det**0.5) ** 2 == icosian_det,
        },
        "interpretation": (
            "With the current trace-polar convention, the icosian Z-lattice is "
            "positive definite and even but not unimodular. Its determinant is "
            "625 and its discriminant group has Smith factors 5,5,5,5. The "
            "G0 trace lattice is the orthogonal double with determinant 5^8."
        ),
    }


def write_json(path: Path, data: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="../../3. Certificates/golden_octonions/g3")
    args = parser.parse_args()
    out = Path(args.out).resolve()
    certificate = build_certificate()
    write_json(out / "e8e8_unimodularity_certificate.json", certificate)
    print(json.dumps({
        "icosian_gram_determinant": certificate["icosian_gram_determinant"],
        "g0_gram_determinant": certificate["g0_gram_determinant"],
        "expected_det_1_confirmed": certificate["decision_gate"]["expected_det_1_confirmed"],
        "normalization_risk_triggered": certificate["decision_gate"]["normalization_risk_triggered"],
    }, indent=2, sort_keys=True))
    print(f"E8/E8 trace-lattice certificate written to {out}")


if __name__ == "__main__":
    main()
