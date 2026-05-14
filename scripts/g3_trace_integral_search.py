from __future__ import annotations

from collections import Counter, defaultdict
import argparse
import json
from pathlib import Path

from golden_octonions import (
    K,
    ZERO,
    ONE,
    HALF,
    PHI,
    Octonion,
    oct_basis_icosian_double,
    h4_quaternion_roots,
    o_coordinates,
)
from g3_search import (
    canonical_line,
    compute_structure,
    scalar_mul_vector,
    vector_inner,
)


def trace_is_integer(x: K) -> bool:
    return x.trace_q().denominator == 1


def trace_text(x: K) -> str:
    return str(x.trace_q())


def basis_vectors(n: int = 8) -> list[list[K]]:
    return [[ONE if i == j else ZERO for i in range(n)] for j in range(n)]


def vector_scale_phi(v: list[K]) -> list[K]:
    return scalar_mul_vector(PHI, v)


def vector_mod2_line(coords: list[K]) -> tuple[int, ...]:
    codes = []
    for c in coords:
        if not c.is_integral():
            raise ValueError(f"non-integral H4 coordinate: {c}")
        codes.append((int(c.a) & 1) | ((int(c.b) & 1) << 1))
    return canonical_line(tuple(codes))


def module_trace_integral_generators(v: list[K]) -> list[list[K]]:
    e = basis_vectors(8)
    return e + [vector_scale_phi(x) for x in e] + [v, vector_scale_phi(v)]


def module_trace_integral(v: list[K], gram) -> tuple[bool, dict | None]:
    gens = module_trace_integral_generators(v)
    labels = [f"e{i}" for i in range(8)] + [f"phi*e{i}" for i in range(8)] + ["v", "phi*v"]
    for i, x in enumerate(gens):
        value = vector_inner(x, x, gram)
        if not trace_is_integer(value):
            return False, {
                "kind": "trace_norm_failure",
                "element": labels[i],
                "norm_value": str(value),
                "trace_norm": trace_text(value),
            }
    for i, x in enumerate(gens):
        for j in range(i + 1, len(gens)):
            y = gens[j]
            value = 2 * vector_inner(x, y, gram)
            if not trace_is_integer(value):
                return False, {
                    "kind": "polar_pairing_failure",
                    "left": labels[i],
                    "right": labels[j],
                    "polar_pairing": str(value),
                    "trace": trace_text(value),
                }
    return True, None


def h4_mixed_half_root_records():
    basis = oct_basis_icosian_double()
    h4 = h4_quaternion_roots()
    zero_q = basis[0].b
    zero_first_q = basis[0].a.scale(0)
    first_coords = [o_coordinates(Octonion(r, zero_q), basis) for r in h4]
    second_coords = [o_coordinates(Octonion(zero_first_q, r), basis) for r in h4]
    for coords_a in first_coords:
        for coords_b in second_coords:
            coords = coords_a[:4] + coords_b[4:]
            yield {
                "coords": coords,
                "line": vector_mod2_line(coords),
                "v": scalar_mul_vector(HALF, coords),
            }


def run_h4_ansatz_search() -> tuple[dict, dict]:
    basis = oct_basis_icosian_double()
    _, _, gram = compute_structure(basis)
    e = basis_vectors(8)
    raw_pair_count = 0
    raw_line_counts: dict[tuple[int, ...], int] = defaultdict(int)
    line_representatives: dict[tuple[int, ...], list[K]] = {}
    raw_counts = Counter()
    projective_counts = Counter()
    first_examples: dict[str, dict] = {}

    for record in h4_mixed_half_root_records():
        raw_pair_count += 1
        line = record["line"]
        raw_line_counts[line] += 1
        line_representatives.setdefault(line, record["coords"])
        v = record["v"]
        norm = vector_inner(v, v, gram)
        self_ok = trace_is_integer(norm)
        g0_polar_ok = all(trace_is_integer(2 * vector_inner(v, x, gram)) for x in e)
        if self_ok:
            raw_counts["passing_self_norm"] += 1
        else:
            first_examples.setdefault("self_norm_failure", {
                "norm": str(norm),
                "trace_norm": trace_text(norm),
                "line": list(line),
            })
        if g0_polar_ok:
            raw_counts["passing_pairings_g0"] += 1

    for line, coords in line_representatives.items():
        v = scalar_mul_vector(HALF, coords)
        norm = vector_inner(v, v, gram)
        if trace_is_integer(norm):
            projective_counts["passing_self_norm"] += 1
        if all(trace_is_integer(2 * vector_inner(v, x, gram)) for x in e):
            projective_counts["passing_pairings_g0"] += 1
        module_ok, failure = module_trace_integral(v, gram)
        if module_ok:
            projective_counts["passing_module_trace_integral"] += 1
        else:
            first_examples.setdefault("module_trace_integrality_failure", {
                "line": list(line),
                "failure": failure,
                "trace_norm_v": trace_text(norm),
                "trace_norm_phi_v": trace_text(vector_inner(vector_scale_phi(v), vector_scale_phi(v), gram)),
                "norm_v": str(norm),
                "norm_phi_v": str(vector_inner(vector_scale_phi(v), vector_scale_phi(v), gram)),
            })

    g0_pairing_lines = [
        line
        for line, coords in line_representatives.items()
        if all(trace_is_integer(2 * vector_inner(scalar_mul_vector(HALF, coords), x, gram)) for x in e)
    ]

    mutual_ok_lines = []
    integral_mutual_pair_count = 0
    total_mutual_pair_count = 0
    first_mutual_failure = None
    for i, line in enumerate(g0_pairing_lines):
        v = scalar_mul_vector(HALF, line_representatives[line])
        line_ok = True
        for j, other_line in enumerate(g0_pairing_lines):
            w = scalar_mul_vector(HALF, line_representatives[other_line])
            value = 2 * vector_inner(v, w, gram)
            ok = trace_is_integer(value)
            if i <= j:
                total_mutual_pair_count += 1
                if ok:
                    integral_mutual_pair_count += 1
            if not ok:
                line_ok = False
                if first_mutual_failure is None:
                    first_mutual_failure = {
                        "line": list(line),
                        "other_line": list(other_line),
                        "polar_pairing": str(value),
                        "trace": trace_text(value),
                    }
        if line_ok:
            mutual_ok_lines.append(line)

    module_survivor_count = projective_counts["passing_module_trace_integral"]
    result = {
        "search": "G3-B-prime H4 mixed half-root ansatz",
        "convention": {
            "candidate": "v = (a+b*l)/2 with a,b in H4",
            "self_trace_norm": "Tr_{K/Q}(N(v)) in Z",
            "roadmap_polar_pairing_filter": "Tr_{K/Q}(B(v,e_i)) in Z for B=2*<,>",
            "module_trace_integrality_check": "Tr_{K/Q}(N(x)) in Z and Tr_{K/Q}(B(x,y)) in Z for the full Z-basis of G0 + Z[phi]*v",
        },
        "raw_pair_count": raw_pair_count,
        "projective_coset_count": len(line_representatives),
        "raw_pairs_per_projective_coset_distribution": dict(sorted(Counter(raw_line_counts.values()).items())),
        "passing_self_norm_count": raw_counts["passing_self_norm"],
        "passing_self_norm_projective_count": projective_counts["passing_self_norm"],
        "passing_pairings_g0_count": raw_counts["passing_pairings_g0"],
        "passing_pairings_g0_projective_count": projective_counts["passing_pairings_g0"],
        "passing_mutual_pairings_count": len(mutual_ok_lines),
        "integral_mutual_pair_count": integral_mutual_pair_count,
        "total_mutual_pair_count": total_mutual_pair_count,
        "passing_module_trace_integral_projective_count": module_survivor_count,
        "passing_module_trace_integral_raw_count": sum(raw_line_counts[line] for line, coords in line_representatives.items() if module_trace_integral(scalar_mul_vector(HALF, coords), gram)[0]),
        "multiplicative_closure_depth": None,
        "multiplicative_closure_status": "not_run_no_trace_integral_Zphi_module_survivors",
        "module_trace_integral_after_closure": False,
        "conjugation_closed_after_closure": None,
        "shell_size_trace_norm_one": 0,
        "mixed_root_count": 0,
        "projective_gram_graph_components": [],
        "cartan_coefficients_in_zphi": None,
        "is_indecomposable": False,
        "first_examples": {
            **first_examples,
            "mutual_pairing_failure": first_mutual_failure,
        },
        "decision": (
            "The H4 mixed half-root ansatz has no Z[phi]-module survivor under "
            "the stated trace-integral convention. Each raw half-root has "
            "Tr(N(v))=1, but adjoining the full Z[phi]-line also adjoins phi*v, "
            "whose trace norm is 3/2 in the first failure example."
        ),
    }
    summary = {
        "raw_pair_count": result["raw_pair_count"],
        "projective_coset_count": result["projective_coset_count"],
        "passing_self_norm_count": result["passing_self_norm_count"],
        "passing_pairings_g0_count": result["passing_pairings_g0_count"],
        "passing_pairings_g0_projective_count": result["passing_pairings_g0_projective_count"],
        "passing_mutual_pairings_count": result["passing_mutual_pairings_count"],
        "passing_module_trace_integral_projective_count": result["passing_module_trace_integral_projective_count"],
        "multiplicative_closure_depth": result["multiplicative_closure_depth"],
        "module_trace_integral_after_closure": result["module_trace_integral_after_closure"],
        "shell_size_trace_norm_one": result["shell_size_trace_norm_one"],
        "mixed_root_count": result["mixed_root_count"],
        "projective_gram_graph_components": result["projective_gram_graph_components"],
        "cartan_coefficients_in_zphi": result["cartan_coefficients_in_zphi"],
        "is_indecomposable": result["is_indecomposable"],
        "decision": result["decision"],
    }
    return result, summary


def write_json(path: Path, data: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="../../3. Certificates/golden_octonions/g3_b_prime")
    args = parser.parse_args()
    out = Path(args.out).resolve()
    result, summary = run_h4_ansatz_search()
    write_json(out / "h4_ansatz_search.json", result)
    write_json(out / "h4_ansatz_summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"G3-B' H4 ansatz certificates written to {out}")


if __name__ == "__main__":
    main()
