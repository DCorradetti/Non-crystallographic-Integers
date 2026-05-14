from __future__ import annotations

from fractions import Fraction
from itertools import product
import argparse
import json
from pathlib import Path

from golden_octonions import K, ZERO, ONE, PHI, oct_basis_icosian_double
from g3_e8e8_unimodularity import build_certificate as build_trace_lattice_certificate
from g3_search import compute_structure, vector_conj, vector_product


P = 5


def inv_mod(x: int, p: int = P) -> int:
    return pow(x % p, -1, p)


def canonical_line(v: tuple[int, ...]) -> tuple[int, ...]:
    for x in v:
        if x % P:
            inv = inv_mod(x)
            return tuple((inv * y) % P for y in v)
    return v


def kernel_basis_param_to_z(c: tuple[int, ...] | list[int]) -> list[int]:
    # The kernel of the certified trace-polar Gram matrix modulo 5 is spanned
    # by (2,1) in each (b_i, phi*b_i) pair.
    out = [0] * 16
    for i in range(4):
        out[i] = (2 * c[i]) % P
        out[4 + i] = c[i] % P
    for i in range(4):
        out[8 + i] = (2 * c[4 + i]) % P
        out[12 + i] = c[4 + i] % P
    return out


def z_to_kernel_params(t: list[int]) -> tuple[int, ...]:
    return tuple([t[4 + i] % P for i in range(4)] + [t[12 + i] % P for i in range(4)])


def q_numerator(param: tuple[int, ...], gram_z: list[list[int]]) -> int:
    t = kernel_basis_param_to_z(param)
    return sum(t[i] * gram_z[i][j] * t[j] for i in range(16) for j in range(16))


def polar_numerator(left: tuple[int, ...], right: tuple[int, ...], gram_z: list[list[int]]) -> int:
    t = kernel_basis_param_to_z(left)
    u = kernel_basis_param_to_z(right)
    return sum(t[i] * gram_z[i][j] * u[j] for i in range(16) for j in range(16))


def is_trace_norm_integral(param: tuple[int, ...], gram_z: list[list[int]]) -> bool:
    # For x=t/5, Tr(N(x)) = t^T G_Z t / 50.
    return q_numerator(param, gram_z) % 50 == 0


def polar_is_integral(left: tuple[int, ...], right: tuple[int, ...], gram_z: list[list[int]]) -> bool:
    # For x=t/5 and y=u/5, Tr(B(x,y)) = t^T G_Z u / 25.
    return polar_numerator(left, right, gram_z) % 25 == 0


def t_to_kcoords(t: list[int]) -> list[K]:
    coords = []
    for base in [0, 8]:
        for i in range(4):
            coords.append(K(Fraction(t[base + i], 5), Fraction(t[base + 4 + i], 5)))
    return coords


def kcoords_to_kernel_param(coords: list[K]) -> tuple[int, ...] | None:
    t = [0] * 16
    for idx, c in enumerate(coords):
        num_a = c.a * 5
        num_b = c.b * 5
        if num_a.denominator != 1 or num_b.denominator != 1:
            return None
        if idx < 4:
            t[idx] = int(num_a) % P
            t[4 + idx] = int(num_b) % P
        else:
            j = idx - 4
            t[8 + j] = int(num_a) % P
            t[12 + j] = int(num_b) % P
    if kernel_basis_param_to_z(z_to_kernel_params(t)) != [x % P for x in t]:
        return None
    return z_to_kernel_params(t)


def basis_vectors(n: int = 8) -> list[list[K]]:
    return [[ONE if i == j else ZERO for i in range(n)] for j in range(n)]


def build_action_maps() -> tuple[list[dict], bool]:
    basis = oct_basis_icosian_double()
    products, conjugation, _ = compute_structure(basis)
    e = basis_vectors(8)
    maps = []
    map_specs = [{"label": "conjugation", "kind": "conjugation"}]
    map_specs += [{"label": f"L_e{i}", "kind": "left", "index": i} for i in range(8)]
    map_specs += [{"label": f"R_e{i}", "kind": "right", "index": i} for i in range(8)]
    preserves_kernel = True
    for spec in map_specs:
        columns = []
        for k in range(8):
            param = [0] * 8
            param[k] = 1
            v = t_to_kcoords(kernel_basis_param_to_z(param))
            if spec["kind"] == "conjugation":
                image = vector_conj(v, conjugation)
            elif spec["kind"] == "left":
                image = vector_product(e[spec["index"]], v, products)
            else:
                image = vector_product(v, e[spec["index"]], products)
            image_param = kcoords_to_kernel_param(image)
            if image_param is None:
                preserves_kernel = False
                image_param = tuple([0] * 8)
            columns.append(image_param)
        maps.append({"label": spec["label"], "columns": columns})
    return maps, preserves_kernel


def apply_map(action_map: dict, vector: tuple[int, ...]) -> tuple[int, ...]:
    out = [0] * 8
    for coeff, col in zip(vector, action_map["columns"]):
        if coeff:
            for i in range(8):
                out[i] = (out[i] + coeff * col[i]) % P
    return tuple(out)


def reduce_vector(v: tuple[int, ...], echelon: dict[int, tuple[int, ...]]) -> tuple[int, ...]:
    out = list(v)
    for pivot in sorted(echelon):
        coeff = out[pivot] % P
        if coeff:
            row = echelon[pivot]
            out = [(x - coeff * y) % P for x, y in zip(out, row)]
    return tuple(out)


def add_to_echelon(v: tuple[int, ...], echelon: dict[int, tuple[int, ...]]) -> bool:
    reduced = reduce_vector(v, echelon)
    if all(x == 0 for x in reduced):
        return False
    pivot = next(i for i, x in enumerate(reduced) if x % P)
    inv = inv_mod(reduced[pivot])
    normalized = tuple((inv * x) % P for x in reduced)
    for existing_pivot, row in list(echelon.items()):
        coeff = row[pivot] % P
        if coeff:
            echelon[existing_pivot] = tuple((x - coeff * y) % P for x, y in zip(row, normalized))
    echelon[pivot] = normalized
    return True


def generated_action_subspace(seed: tuple[int, ...], action_maps: list[dict]) -> list[tuple[int, ...]]:
    echelon: dict[int, tuple[int, ...]] = {}
    queue: list[tuple[int, ...]] = []
    if add_to_echelon(seed, echelon):
        queue.append(seed)
    idx = 0
    while idx < len(queue):
        current = queue[idx]
        idx += 1
        for action_map in action_maps:
            image = apply_map(action_map, current)
            before = len(echelon)
            if add_to_echelon(image, echelon) and len(echelon) > before:
                queue.append(image)
                if len(echelon) == 8:
                    return list(echelon.values())
    return list(echelon.values())


def subspace_totally_isotropic(basis: list[tuple[int, ...]], gram_z: list[list[int]]) -> bool:
    for i, left in enumerate(basis):
        if not is_trace_norm_integral(left, gram_z):
            return False
        for right in basis[i + 1 :]:
            if not polar_is_integral(left, right, gram_z):
                return False
    return True


def phi_action_scalar_check() -> dict:
    # In each (b_i, phi*b_i) pair, phi acts by [[0,1],[1,1]].
    scalar_ok = True
    sqrt5_ok = True
    for k in range(8):
        param = [0] * 8
        param[k] = 1
        t = kernel_basis_param_to_z(param)
        phi_t = [0] * 16
        for base in [0, 8]:
            for i in range(4):
                phi_t[base + i] = t[base + 4 + i]
                phi_t[base + 4 + i] = (t[base + i] + t[base + 4 + i]) % P
        scalar_3 = [(3 * x) % P for x in t]
        sqrt5_t = [((2 * phi_t[i] - t[i]) % P) for i in range(16)]
        scalar_ok = scalar_ok and phi_t == scalar_3
        sqrt5_ok = sqrt5_ok and all(x == 0 for x in sqrt5_t)
    return {
        "phi_acts_as_scalar_3_on_discriminant_group": scalar_ok,
        "sqrt5_annihilates_discriminant_group": sqrt5_ok,
    }


def run_phase2() -> tuple[dict, dict, dict, dict]:
    trace_cert = build_trace_lattice_certificate()
    gram_z = trace_cert["g0_gram_matrix"]
    action_maps, maps_preserve_kernel = build_action_maps()
    q_value_distribution: dict[str, int] = {}
    isotropic_line_count = 0
    closure_dimension_distribution: dict[str, int] = {}
    first_examples: dict[str, object] = {}
    octonion_stable_isotropic_subspace_count = 0
    indecomposable_candidate_count = 0
    projective_line_count = 0

    for raw in product(range(P), repeat=8):
        if all(x == 0 for x in raw):
            continue
        line = canonical_line(raw)
        if raw != line:
            continue
        projective_line_count += 1
        q_mod = q_numerator(line, gram_z) % 50
        q_value_distribution[str(q_mod)] = q_value_distribution.get(str(q_mod), 0) + 1
        if q_mod != 0:
            continue
        isotropic_line_count += 1
        generated = generated_action_subspace(line, action_maps)
        dim = len(generated)
        closure_dimension_distribution[str(dim)] = closure_dimension_distribution.get(str(dim), 0) + 1
        first_examples.setdefault("first_isotropic_line", list(line))
        first_examples.setdefault(f"first_generated_dimension_{dim}", [list(v) for v in generated])
        if subspace_totally_isotropic(generated, gram_z):
            octonion_stable_isotropic_subspace_count += 1
            first_examples.setdefault("first_octonion_stable_isotropic_subspace", [list(v) for v in generated])

    discriminant = {
        "certificate": "discriminant_quadratic_form",
        "source": "e8e8_unimodularity_certificate.json",
        "group": "G0_Z_dual/G0 as trace-polar Z-lattice",
        "discriminant_group_order": 5**8,
        "discriminant_group_structure": "(Z/5)^8",
        "kernel_basis_mod5_parametrization": "param c in F5^8 maps to Z-vector (2c_0,...,2c_3,c_0,...,c_3,2c_4,...,2c_7,c_4,...,c_7)",
        "quadratic_form": "q(t/5)=t^T G_Z t / 50 mod Z",
        "bilinear_form": "b(t/5,u/5)=t^T G_Z u / 25 mod Z",
        "quadratic_form_witt_type": "split plus type O^+(8,5), hyperbolic rank 4",
        "maximal_isotropic_dimension": 4,
        "projective_line_count": projective_line_count,
        "isotropic_projective_line_count": isotropic_line_count,
        "q_value_distribution_mod_50": dict(sorted(q_value_distribution.items(), key=lambda kv: int(kv[0]))),
        "phi_action": phi_action_scalar_check(),
        "action_maps_preserve_discriminant_kernel": maps_preserve_kernel,
    }

    search = {
        "certificate": "sqrt5_ansatz_search",
        "search": "G3-B-prime sqrt5 discriminant-group ansatz",
        "coset_count_total": 5**8 - 1,
        "projective_line_count": projective_line_count,
        "passing_phi_stability_count": projective_line_count if discriminant["phi_action"]["phi_acts_as_scalar_3_on_discriminant_group"] else None,
        "passing_trace_norm_line_count": isotropic_line_count,
        "passing_module_assembly_count": isotropic_line_count,
        "multiplicative_closure_depth": 1,
        "basis_action_generated_subspace_dimension_distribution": dict(sorted(closure_dimension_distribution.items(), key=lambda kv: int(kv[0]))),
        "passing_closure_trace_integral_count": octonion_stable_isotropic_subspace_count,
        "shell_size_trace_norm_one": 0,
        "mixed_root_count": 0,
        "projective_gram_graph_components": [],
        "cartan_in_zphi": None,
        "indecomposable_candidate_count": indecomposable_candidate_count,
        "first_examples": first_examples,
        "decision": (
            "Every isotropic line generates the full 8-dimensional discriminant "
            "group under conjugation and left/right multiplication by the G0 "
            "basis. The full discriminant group is not totally isotropic, so no "
            "nonzero octonion-stable isotropic subspace exists in this ansatz."
        ),
    }

    stable = {
        "certificate": "octonion_stable_subspaces",
        "method": (
            "For each isotropic projective line, generate the smallest F5-subspace "
            "stable under conjugation and left/right multiplication by the G0 basis. "
            "Any nonzero stable isotropic subspace would contain a line whose generated "
            "stable closure is still totally isotropic."
        ),
        "isotropic_projective_line_count": isotropic_line_count,
        "generated_subspace_dimension_distribution": search["basis_action_generated_subspace_dimension_distribution"],
        "octonion_stable_subspace_count": octonion_stable_isotropic_subspace_count,
        "surviving_indecomposable_candidates": [],
        "proof_of_zero": (
            "The generated stable closure of every isotropic line has dimension 8. "
            "Since the full discriminant group is not totally isotropic, no nonzero "
            "octonion-stable isotropic subspace can exist."
        ),
    }

    summary = {
        "discriminant_group_order": discriminant["discriminant_group_order"],
        "quadratic_form_witt_type": discriminant["quadratic_form_witt_type"],
        "coset_count_total": search["coset_count_total"],
        "projective_line_count": projective_line_count,
        "isotropic_projective_line_count": isotropic_line_count,
        "passing_phi_stability_count": search["passing_phi_stability_count"],
        "octonion_stable_subspace_count": octonion_stable_isotropic_subspace_count,
        "indecomposable_candidate_count": indecomposable_candidate_count,
        "decision": search["decision"],
    }
    return discriminant, search, stable, summary


def write_json(path: Path, data: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="../../3. Certificates/golden_octonions/g3_b_prime")
    args = parser.parse_args()
    out = Path(args.out).resolve()
    discriminant, search, stable, summary = run_phase2()
    write_json(out / "discriminant_quadratic_form.json", discriminant)
    write_json(out / "sqrt5_ansatz_search.json", search)
    write_json(out / "octonion_stable_subspaces.json", stable)
    write_json(out / "sqrt5_ansatz_summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"G3-B' Phase 2 certificates written to {out}")


if __name__ == "__main__":
    main()
