from __future__ import annotations

import argparse
import json
from pathlib import Path

from golden_octonions import K, ZERO, ONE, TWO, oct_basis_icosian_double
from g3_search import compute_structure


BASIS_LABELS = ["1", "i", "h", "g", "l", "il", "hl", "gl"]


def k_payload(x: K) -> dict:
    return {
        "text": str(x),
        "coefficients_in_1_phi": x.as_pair(),
        "field_norm_to_Q": str(x.norm_q()),
        "is_integral": x.is_integral(),
        "is_unit": x.is_integral() and abs(x.norm_q()) == 1,
    }


def k_text_matrix(matrix: list[list[K]]) -> list[list[str]]:
    return [[str(x) for x in row] for row in matrix]


def determinant(matrix: list[list[K]]) -> K:
    n = len(matrix)
    a = [row[:] for row in matrix]
    det = ONE
    sign = 1
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if not a[row][col].is_zero():
                pivot = row
                break
        if pivot is None:
            return ZERO
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            sign *= -1
        pivot_value = a[col][col]
        det = det * pivot_value
        pivot_inv = pivot_value.inverse()
        for row in range(col + 1, n):
            factor = a[row][col] * pivot_inv
            if factor.is_zero():
                continue
            for idx in range(col, n):
                a[row][idx] = a[row][idx] - factor * a[col][idx]
    return det if sign == 1 else -det


def inverse_matrix(matrix: list[list[K]]) -> list[list[K]]:
    n = len(matrix)
    a = [
        row[:] + [ONE if i == j else ZERO for j in range(n)]
        for i, row in enumerate(matrix)
    ]
    pivot_row = 0
    for col in range(n):
        pivot = None
        for row in range(pivot_row, n):
            if not a[row][col].is_zero():
                pivot = row
                break
        if pivot is None:
            raise ValueError("matrix is singular")
        a[pivot_row], a[pivot] = a[pivot], a[pivot_row]
        pivot_inv = a[pivot_row][col].inverse()
        a[pivot_row] = [pivot_inv * x for x in a[pivot_row]]
        for row in range(n):
            if row == pivot_row:
                continue
            factor = a[row][col]
            if factor.is_zero():
                continue
            a[row] = [x - factor * y for x, y in zip(a[row], a[pivot_row])]
        pivot_row += 1
    return [row[n:] for row in a]


def matrix_product(a: list[list[K]], b: list[list[K]]) -> list[list[K]]:
    rows = len(a)
    cols = len(b[0])
    mid = len(b)
    out = [[ZERO for _ in range(cols)] for _ in range(rows)]
    for i in range(rows):
        for j in range(cols):
            total = ZERO
            for k in range(mid):
                total = total + a[i][k] * b[k][j]
            out[i][j] = total
    return out


def identity_matrix(n: int) -> list[list[K]]:
    return [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]


def matrix_equal(a: list[list[K]], b: list[list[K]]) -> bool:
    return all(x == y for row_a, row_b in zip(a, b) for x, y in zip(row_a, row_b))


def matrix_symmetric(matrix: list[list[K]]) -> bool:
    return all(matrix[i][j] == matrix[j][i] for i in range(len(matrix)) for j in range(len(matrix)))


def build_certificate() -> dict:
    basis = oct_basis_icosian_double()
    _, _, inner_matrix = compute_structure(basis)
    polar_matrix = [[TWO * x for x in row] for row in inner_matrix]
    inner_det = determinant(inner_matrix)
    polar_det = determinant(polar_matrix)
    polar_inverse = inverse_matrix(polar_matrix)
    identity = identity_matrix(len(polar_matrix))
    inverse_left_ok = matrix_equal(matrix_product(polar_inverse, polar_matrix), identity)
    inverse_right_ok = matrix_equal(matrix_product(polar_matrix, polar_inverse), identity)
    polar_inverse_integral = all(x.is_integral() for row in polar_inverse for x in row)
    polar_matrix_integral = all(x.is_integral() for row in polar_matrix for x in row)
    unimodular = polar_det.is_integral() and abs(polar_det.norm_q()) == 1 and polar_inverse_integral
    return {
        "certificate": "dual_discriminant_certificate",
        "order": "G0 = I + I*l",
        "coefficient_ring": "R = Z[phi], phi^2 = phi + 1",
        "basis_labels": BASIS_LABELS,
        "pairing_convention": {
            "inner_product": "<x,y> = Re(x*conj(y))",
            "norm_polar_bilinear": "B(x,y) = N(x+y)-N(x)-N(y) = 2*<x,y>",
            "dual_used_for_strict_norm_integral_overorders": "B",
        },
        "inner_product_matrix": k_text_matrix(inner_matrix),
        "inner_product_determinant": k_payload(inner_det),
        "norm_polar_pairing_matrix": k_text_matrix(polar_matrix),
        "norm_polar_determinant": k_payload(polar_det),
        "norm_polar_inverse_matrix": k_text_matrix(polar_inverse),
        "dual_basis_coordinates_in_G0_basis": k_text_matrix(polar_inverse),
        "ideal_factorization": {
            "norm_polar_determinant": "unit ideal",
            "candidate_prime_ideals": [],
            "reason": "The norm-polar determinant is a unit of Z[phi].",
        },
        "quotient_invariants": {
            "G0_dual_over_G0": "trivial",
            "cardinality": 1,
            "local_components": [],
        },
        "checks": {
            "inner_product_matrix_symmetric": matrix_symmetric(inner_matrix),
            "norm_polar_matrix_symmetric": matrix_symmetric(polar_matrix),
            "norm_polar_matrix_entries_integral": polar_matrix_integral,
            "norm_polar_inverse_entries_integral": polar_inverse_integral,
            "inverse_left_verified": inverse_left_ok,
            "inverse_right_verified": inverse_right_ok,
            "determinant_verified": polar_det == K(2, 3),
            "norm_polar_determinant_is_unit": abs(polar_det.norm_q()) == 1,
            "G0_is_unimodular_for_norm_polar_pairing": unimodular,
        },
        "conclusion": (
            "The norm-polar dual G0# equals G0. Therefore a strict overorder "
            "of G0 with integral norms cannot be obtained by adjoining any "
            "fractional ideal layer inside K tensor_R G0."
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
    write_json(out / "dual_discriminant_certificate.json", certificate)
    print(json.dumps({
        "norm_polar_determinant": certificate["norm_polar_determinant"]["text"],
        "norm_polar_determinant_norm": certificate["norm_polar_determinant"]["field_norm_to_Q"],
        "G0_is_unimodular_for_norm_polar_pairing": certificate["checks"][
            "G0_is_unimodular_for_norm_polar_pairing"
        ],
        "candidate_prime_ideals": certificate["ideal_factorization"]["candidate_prime_ideals"],
    }, indent=2, sort_keys=True))
    print(f"Dual/discriminant certificate written to {out}")


if __name__ == "__main__":
    main()
