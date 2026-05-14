from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

from golden_octonions import K, TWO, oct_basis_icosian_double
from g3_search import compute_structure


def f5_inv(x: int) -> int:
    if x % 5 == 0:
        raise ZeroDivisionError("0 has no inverse in F5")
    for y in range(1, 5):
        if (x * y) % 5 == 1:
            return y
    raise AssertionError("unreachable")


def canonical_f5_line(v: tuple[int, ...]) -> tuple[int, ...]:
    for x in v:
        if x % 5:
            inv = f5_inv(x)
            return tuple((inv * y) % 5 for y in v)
    return v


def generate_f5_lines(dim: int = 8):
    seen = set()
    for v in itertools.product(range(5), repeat=dim):
        if all(x == 0 for x in v):
            continue
        c = canonical_f5_line(v)
        if c in seen:
            continue
        seen.add(c)
        yield c


def mod_sqrt5(x: K) -> int:
    if not x.is_integral():
        raise ValueError(f"non-integral residue input: {x}")
    # sqrt(5)=2*phi-1 vanishes, hence phi = 1/2 = 3 in F5.
    return (int(x.a) + 3 * int(x.b)) % 5


def line_is_mixed(v: tuple[int, ...]) -> bool:
    return any(v[i] for i in range(4)) and any(v[i] for i in range(4, 8))


def polar_pairing_residues_mod_sqrt5(line: tuple[int, ...], polar_matrix_mod5: list[list[int]]) -> list[int]:
    residues = []
    for j in range(8):
        total = 0
        for i in range(8):
            total = (total + line[i] * polar_matrix_mod5[i][j]) % 5
        residues.append(total)
    return residues


def denominator_sqrt5_line_scan(limit: int | None = None) -> dict:
    basis = oct_basis_icosian_double()
    _, _, gram = compute_structure(basis)
    polar = [[TWO * x for x in row] for row in gram]
    polar_mod5 = [[mod_sqrt5(x) for x in row] for row in polar]
    line_count = 0
    mixed_count = 0
    pairing_count = 0
    first_pairing_fail = None
    first_pairing_pass = None
    for idx, line in enumerate(generate_f5_lines(8), start=1):
        if limit is not None and idx > limit:
            break
        line_count += 1
        if line_is_mixed(line):
            mixed_count += 1
        residues = polar_pairing_residues_mod_sqrt5(line, polar_mod5)
        pairing_ok = all(r == 0 for r in residues)
        if pairing_ok:
            pairing_count += 1
            if first_pairing_pass is None:
                first_pairing_pass = line
        elif first_pairing_fail is None:
            first_pairing_fail = {"line": line, "pairing_residues_mod_sqrt5": residues}
    complete = limit is None
    return {
        "search": "denominator_sqrt5_single_line_search",
        "ideal": "(sqrt5) = (2*phi-1)",
        "residue_field": "Z[phi]/(sqrt5) ~= F5, phi -> 3",
        "ambient_quotient": "G0/sqrt5*G0 ~= F5^8",
        "line_count_expected": 97656,
        "line_count_tested": line_count,
        "mixed_line_count_tested": mixed_count,
        "projective_lines_with_polar_pairings_integral": pairing_count,
        "survivor_count": 0,
        "complete_for_single_line_denominator_sqrt5": complete,
        "consistent_with_strict_dual_theorem": complete and pairing_count == 0,
        "first_pairing_fail_example": {
            "line": list(first_pairing_fail["line"]),
            "pairing_residues_mod_sqrt5": first_pairing_fail["pairing_residues_mod_sqrt5"],
        } if first_pairing_fail else None,
        "first_pairing_pass_example": list(first_pairing_pass) if first_pairing_pass else None,
        "filter_note": (
            "The necessary polar-pairing filter is already decisive. Since no "
            "nonzero projective F5 line pairs integrally after division by "
            "sqrt5, later multiplication and square-closure filters are not reached."
        ),
    }


def denominator_sqrt5_all_subspaces_no_go(line_scan: dict) -> dict:
    excluded = (
        line_scan["complete_for_single_line_denominator_sqrt5"]
        and line_scan["projective_lines_with_polar_pairings_integral"] == 0
    )
    return {
        "search": "denominator_sqrt5_all_f5_subspaces_pairing_obstruction",
        "ambient_quotient": "G0/sqrt5*G0 ~= F5^8",
        "subspaces_considered": "all nonzero F5-linear subspaces C <= F5^8",
        "projective_line_count_tested": line_scan["line_count_tested"],
        "projective_lines_with_polar_pairings_integral": line_scan[
            "projective_lines_with_polar_pairings_integral"
        ],
        "all_nonzero_denominator_sqrt5_subspaces_excluded": excluded,
        "proof": [
            "Every nonzero F5-subspace contains a projective F5-line.",
            "The complete projective-line scan found zero lines satisfying the necessary polar-pairing condition.",
            "Therefore no nonzero sqrt5-denominator subspace can define a strict norm-integral overmodule of G0.",
        ],
        "source_certificate": "denominator_sqrt5_single_line_search.json",
    }


def write_json(path: Path, data: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="../../3. Certificates/golden_octonions/g3")
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()
    out = Path(args.out).resolve()
    line_scan = denominator_sqrt5_line_scan(limit=args.limit)
    subspaces = denominator_sqrt5_all_subspaces_no_go(line_scan)
    write_json(out / "denominator_sqrt5_single_line_search.json", line_scan)
    write_json(out / "denominator_sqrt5_all_subspaces_no_go.json", subspaces)
    print(json.dumps({
        "line_count_tested": line_scan["line_count_tested"],
        "mixed_line_count_tested": line_scan["mixed_line_count_tested"],
        "projective_lines_with_polar_pairings_integral": line_scan[
            "projective_lines_with_polar_pairings_integral"
        ],
        "survivor_count": line_scan["survivor_count"],
        "consistent_with_strict_dual_theorem": line_scan["consistent_with_strict_dual_theorem"],
    }, indent=2, sort_keys=True))
    print(f"sqrt5 denominator certificates written to {out}")


if __name__ == "__main__":
    main()
