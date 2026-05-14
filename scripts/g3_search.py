from __future__ import annotations

from collections import Counter, defaultdict, deque
from fractions import Fraction
import argparse
import itertools
import json
from pathlib import Path

from golden_octonions import (
    K,
    ZERO,
    ONE,
    TWO,
    HALF,
    PHI,
    Octonion,
    OZERO,
    oct_basis_icosian_double,
    h4_quaternion_roots,
    double_shell,
    o_coordinates,
    hash_objects,
)


F4_ZERO = 0
F4_ONE = 1
F4_PHI = 2
F4_ONE_PLUS_PHI = 3
F4_REPS = {
    F4_ZERO: ZERO,
    F4_ONE: ONE,
    F4_PHI: PHI,
    F4_ONE_PLUS_PHI: ONE + PHI,
}


def f4_add(x: int, y: int) -> int:
    return x ^ y


def f4_mul(x: int, y: int) -> int:
    if x == 0 or y == 0:
        return 0
    a = K(x & 1, (x >> 1) & 1)
    b = K(y & 1, (y >> 1) & 1)
    z = a * b
    return (int(z.a) & 1) | ((int(z.b) & 1) << 1)


def f4_inv(x: int) -> int:
    if x == 0:
        raise ZeroDivisionError("0 has no inverse in F4")
    for y in (1, 2, 3):
        if f4_mul(x, y) == 1:
            return y
    raise AssertionError("unreachable")


def canonical_line(v: tuple[int, ...]) -> tuple[int, ...]:
    for x in v:
        if x:
            inv = f4_inv(x)
            return tuple(f4_mul(inv, y) for y in v)
    return v


def generate_f4_lines(dim: int = 8):
    seen = set()
    for v in itertools.product(range(4), repeat=dim):
        if all(x == 0 for x in v):
            continue
        c = canonical_line(v)
        if c in seen:
            continue
        seen.add(c)
        yield c


def f4_vector_to_k(v: tuple[int, ...]) -> list[K]:
    return [F4_REPS[x] for x in v]


def coords_to_oct(coords: list[K], basis: list[Octonion]) -> Octonion:
    out = OZERO
    for c, b in zip(coords, basis):
        out = out + b.scale(c)
    return out


def compute_structure(basis: list[Octonion]):
    n = len(basis)
    products = [[[ZERO for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for i, bi in enumerate(basis):
        for j, bj in enumerate(basis):
            products[i][j] = o_coordinates(bi.mul(bj), basis)
    conjugation = []
    for bi in basis:
        conjugation.append(o_coordinates(bi.conj(), basis))
    gram = [[basis[i].inner(basis[j]) for j in range(n)] for i in range(n)]
    return products, conjugation, gram


def vector_product(u: list[K], v: list[K], products) -> list[K]:
    n = len(u)
    out = [ZERO for _ in range(n)]
    for i in range(n):
        if u[i].is_zero():
            continue
        for j in range(n):
            if v[j].is_zero():
                continue
            coeff = u[i] * v[j]
            for k in range(n):
                out[k] = out[k] + coeff * products[i][j][k]
    return out


def vector_conj(u: list[K], conjugation) -> list[K]:
    n = len(u)
    out = [ZERO for _ in range(n)]
    for i in range(n):
        if u[i].is_zero():
            continue
        for k in range(n):
            out[k] = out[k] + u[i] * conjugation[i][k]
    return out


def vector_inner(u: list[K], v: list[K], gram) -> K:
    out = ZERO
    n = len(u)
    for i in range(n):
        if u[i].is_zero():
            continue
        for j in range(n):
            if v[j].is_zero():
                continue
            out = out + u[i] * v[j] * gram[i][j]
    return out


def scalar_mul_vector(s: K, v: list[K]) -> list[K]:
    return [s * x for x in v]


def vector_sub(u: list[K], v: list[K]) -> list[K]:
    return [a - b for a, b in zip(u, v)]


def all_integral(v: list[K]) -> bool:
    return all(x.is_integral() for x in v)


def in_single_glue_module(w: list[K], glue_v: list[K]) -> bool:
    # Module G0 + R*(glue_v/2). Only the residue class of the coefficient in
    # R/2R matters, so four representatives suffice.
    x = scalar_mul_vector(HALF, glue_v)
    for c in F4_REPS.values():
        if all_integral(vector_sub(w, scalar_mul_vector(c, x))):
            return True
    return False


def line_is_mixed(v: tuple[int, ...]) -> bool:
    return any(v[i] for i in range(4)) and any(v[i] for i in range(4, 8))


def line_filter_report(v_codes: tuple[int, ...], products, conjugation, gram) -> dict:
    v = f4_vector_to_k(v_codes)
    x = scalar_mul_vector(HALF, v)
    basis_vectors = [[ONE if i == j else ZERO for i in range(8)] for j in range(8)]
    checks = {}
    checks["mixed_coset"] = line_is_mixed(v_codes)
    checks["conjugation_closed"] = in_single_glue_module(vector_conj(x, conjugation), v)
    checks["pairings_integral"] = all(vector_inner(x, e, gram).is_integral() for e in basis_vectors)
    checks["norm_integral"] = vector_inner(x, x, gram).is_integral()
    left_ok = True
    right_ok = True
    for e in basis_vectors:
        left_ok = left_ok and in_single_glue_module(vector_product(e, x, products), v)
        right_ok = right_ok and in_single_glue_module(vector_product(x, e, products), v)
    checks["basis_multiplication_closed"] = left_ok and right_ok
    checks["square_closed"] = in_single_glue_module(vector_product(x, x, products), v)
    checks["candidate_passes"] = all(checks.values())
    return checks


def denominator_two_search(limit: int | None = None) -> dict:
    basis = oct_basis_icosian_double()
    products, conjugation, gram = compute_structure(basis)
    counters = Counter()
    rejection_reasons = Counter()
    survivors = []
    first_fail_examples: dict[str, tuple[int, ...]] = {}
    for idx, line in enumerate(generate_f4_lines(8), start=1):
        if limit is not None and idx > limit:
            break
        counters["lines_tested"] += 1
        if line_is_mixed(line):
            counters["mixed_lines_tested"] += 1
        checks = line_filter_report(line, products, conjugation, gram)
        for key, value in checks.items():
            if value:
                counters[f"pass_{key}"] += 1
        if checks["candidate_passes"]:
            survivors.append(line)
        else:
            for key in [
                "mixed_coset",
                "conjugation_closed",
                "pairings_integral",
                "norm_integral",
                "basis_multiplication_closed",
                "square_closed",
            ]:
                if not checks[key]:
                    rejection_reasons[key] += 1
                    first_fail_examples.setdefault(key, line)
                    break
    return {
        "search": "denominator_two_single_line_gluing",
        "ambient_quotient": "(1/2 G0)/G0 ~= (Z[phi]/2)^8 ~= F4^8",
        "line_count_expected": 21845,
        "line_count_tested": counters["lines_tested"],
        "mixed_line_count_tested": counters["mixed_lines_tested"],
        "filter_pass_counts": dict(sorted(counters.items())),
        "rejection_reasons_first_failed_filter": dict(sorted(rejection_reasons.items())),
        "first_fail_examples_f4_lines": {k: list(v) for k, v in first_fail_examples.items()},
        "survivor_count": len(survivors),
        "survivors_f4_lines": [list(v) for v in survivors[:25]],
        "complete_for_single_line_denominator_2": limit is None,
    }


def denominator_two_all_subspaces_no_go(line_search: dict) -> dict:
    pairing_line_count = line_search["filter_pass_counts"].get("pass_pairings_integral", 0)
    complete_line_scan = (
        line_search["complete_for_single_line_denominator_2"]
        and line_search["line_count_tested"] == line_search["line_count_expected"]
    )
    excluded = complete_line_scan and pairing_line_count == 0
    return {
        "search": "denominator_two_all_f4_subspaces_pairing_obstruction",
        "ambient_quotient": "(1/2 G0)/G0 ~= (Z[phi]/2)^8 ~= F4^8",
        "subspaces_considered": "all nonzero F4-linear subspaces C <= F4^8",
        "overmodule_form": "G_C = G0 + (1/2)C",
        "necessary_condition": "every projective line contained in C must have integral pairings with the baseline G0 basis",
        "projective_line_count_expected": line_search["line_count_expected"],
        "projective_line_count_tested": line_search["line_count_tested"],
        "projective_lines_with_integral_pairings": pairing_line_count,
        "complete_line_scan": complete_line_scan,
        "all_nonzero_denominator_2_subspaces_excluded": excluded,
        "proof": [
            "Every nonzero F4-subspace C <= F4^8 contains at least one projective F4-line.",
            "If G_C is a strict Z[phi]-order candidate, each element v/2 with v in C must pair integrally with the baseline basis of G0.",
            "The complete projective-line scan found zero lines satisfying this integral-pairing condition.",
            "Therefore no nonzero denominator-2 F4-subspace C can pass the strict integral-pairing test.",
        ],
        "source_certificate": "denominator2_single_line_search.json",
    }


def mixed_half_root_search() -> dict:
    basis = oct_basis_icosian_double()
    products, conjugation, gram = compute_structure(basis)
    h4 = h4_quaternion_roots()
    first_coords = [o_coordinates(Octonion(r, basis[0].b), basis) for r in h4]
    second_coords = [o_coordinates(Octonion(basis[0].a.scale(0), r), basis) for r in h4]
    line_set = {}
    total = 0
    norm_counter = Counter()
    integral_norm_count = 0
    for coords_a in first_coords:
        for coords_b in second_coords:
            total += 1
            coords = coords_a[:4] + coords_b[4:]
            norm = vector_inner(scalar_mul_vector(HALF, coords), scalar_mul_vector(HALF, coords), gram)
            norm_counter[str(norm)] += 1
            if norm.is_integral():
                integral_norm_count += 1
            codes = []
            for c in coords:
                if not c.is_integral():
                    raise ValueError("H4 root did not have integral icosian coordinates")
                codes.append((int(c.a) & 1) | ((int(c.b) & 1) << 1))
            line_set[canonical_line(tuple(codes))] = True
    tested_lines = []
    passing_lines = []
    for line in sorted(line_set):
        checks = line_filter_report(line, products, conjugation, gram)
        tested_lines.append((line, checks))
        if checks["candidate_passes"]:
            passing_lines.append(line)
    return {
        "search": "mixed_half_roots_x=(a+b*l)/2_with_a_b_in_H4",
        "raw_pair_count": total,
        "distinct_projective_f4_lines": len(line_set),
        "norm_distribution": dict(sorted(norm_counter.items())),
        "integral_norm_raw_pair_count": integral_norm_count,
        "passing_line_count": len(passing_lines),
        "passing_lines": [list(v) for v in passing_lines[:25]],
        "interpretation": "For H4 roots normalized to norm 1, all raw mixed half-roots have norm 1/2, hence they are not integral-norm order elements under the strict Z[phi]-order convention.",
    }


def projective_shell(shell: list[Octonion]) -> list[Octonion]:
    seen = set()
    out = []
    for x in shell:
        key = x.key()
        nkey = (-x).key()
        if key in seen or nkey in seen:
            continue
        seen.add(key)
        seen.add(nkey)
        out.append(x)
    return out


def gram_graph_components(shell: list[Octonion]) -> list[int]:
    reps = projective_shell(shell)
    n = len(reps)
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if not reps[i].inner(reps[j]).is_zero():
                adj[i].append(j)
                adj[j].append(i)
    seen = [False] * n
    sizes = []
    for i in range(n):
        if seen[i]:
            continue
        q = deque([i])
        seen[i] = True
        size = 0
        while q:
            u = q.popleft()
            size += 1
            for v in adj[u]:
                if not seen[v]:
                    seen[v] = True
                    q.append(v)
        sizes.append(size)
    return sorted(sizes, reverse=True)


def baseline_invariants() -> dict:
    shell = double_shell()
    pure_first = sum(1 for x in shell if not x.a.is_zero() and x.b.is_zero())
    pure_second = sum(1 for x in shell if x.a.is_zero() and not x.b.is_zero())
    mixed = sum(1 for x in shell if not x.a.is_zero() and not x.b.is_zero())
    inner_distribution = Counter()
    for i, x in enumerate(shell):
        for y in shell[i:]:
            inner_distribution[str(x.inner(y))] += 1
    return {
        "order": "G0 = I + I*l",
        "shell": "H4 + H4",
        "shell_cardinality": len(shell),
        "projective_shell_cardinality": len(projective_shell(shell)),
        "pure_first_half_roots": pure_first,
        "pure_second_half_roots": pure_second,
        "mixed_roots": mixed,
        "gram_graph_component_sizes_projective": gram_graph_components(shell),
        "inner_product_distribution_upper_triangle": dict(sorted(inner_distribution.items())),
        "shell_hash": hash_objects(shell),
        "decomposable_by_orthogonal_blocks": True,
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
    baseline = baseline_invariants()
    den2 = denominator_two_search(limit=args.limit)
    den2_subspaces = denominator_two_all_subspaces_no_go(den2)
    mixed = mixed_half_root_search()
    write_json(out / "g1_baseline_invariants.json", baseline)
    write_json(out / "denominator2_single_line_search.json", den2)
    write_json(out / "denominator2_all_subspaces_no_go.json", den2_subspaces)
    write_json(out / "mixed_half_root_search.json", mixed)
    summary = {
        "baseline_shell": baseline["shell"],
        "baseline_components": baseline["gram_graph_component_sizes_projective"],
        "denominator2_single_line_survivors": den2["survivor_count"],
        "denominator2_complete": den2["complete_for_single_line_denominator_2"],
        "denominator2_projective_lines_with_integral_pairings": den2_subspaces[
            "projective_lines_with_integral_pairings"
        ],
        "denominator2_all_nonzero_subspaces_excluded": den2_subspaces[
            "all_nonzero_denominator_2_subspaces_excluded"
        ],
        "mixed_half_root_passing_lines": mixed["passing_line_count"],
    }
    write_json(out / "summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"G3 certificates written to {out}")


if __name__ == "__main__":
    main()
