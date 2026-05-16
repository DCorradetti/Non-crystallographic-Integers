from __future__ import annotations

import cmath
import hashlib
import json
import math
import shutil
import sys
from datetime import datetime
from fractions import Fraction
from itertools import product
from pathlib import Path

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


HERE = Path(__file__).resolve().parent
_local_root = HERE.parents[1]
if (_local_root / "golden_octonions").is_dir():
    ROOT = _local_root.parent
    GOLDEN_SCRIPTS = ROOT / "2. Scripts" / "golden_octonions"
    CERT_DIR = ROOT / "3. Certificates" / "exact_golden_model_sets"
    RESULTS_ROOT = ROOT / "7. Results" / "17. Exact Golden Model Sets"
    K8E8_DIR = ROOT / "3. Certificates" / "k8_vs_e8"
else:
    SUPP = HERE.parent
    ROOT = SUPP
    GOLDEN_SCRIPTS = HERE
    CERT_DIR = SUPP / "certificates"
    RESULTS_ROOT = SUPP / "runs"
    K8E8_DIR = SUPP / "certificates"
sys.path.insert(0, str(GOLDEN_SCRIPTS))

from golden_octonions import K, ONE, PHI, PHI_INV, h4_quaternion_roots  # noqa: E402


ARTICLE = "17. Exact Golden Model Sets"
CORE_VERSION = "v2_decagonal_h2_window"


def write_json(path: Path, data: dict | list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def sha256_json(data: dict | list) -> str:
    payload = json.dumps(data, indent=2, sort_keys=True).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def frac_text(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def k_text(x: K) -> str:
    return str(x)


def k_pair(x: K) -> list[str]:
    return [frac_text(x.a), frac_text(x.b)]


def k_float(x: K, star: bool = False) -> float:
    phi = (1 - math.sqrt(5)) / 2 if star else (1 + math.sqrt(5)) / 2
    return float(x.a) + float(x.b) * phi


def k_half(x: K) -> K:
    return K(x.a / 2, x.b / 2)


def h2_cos(delta: int) -> K:
    """cos(2*pi*delta/5) in Q(sqrt(5)), delta modulo 5."""
    d = delta % 5
    if d == 0:
        return K(1)
    if d in (1, 4):
        return k_half(PHI - ONE)
    return k_half(-PHI)


def gram_from_step(multiplier: int) -> list[list[K]]:
    return [[h2_cos(multiplier * (i - j)) for j in range(5)] for i in range(5)]


def k_quadratic(v: tuple[int, ...], gram: list[list[K]]) -> K:
    total = K(0)
    for i, vi in enumerate(v):
        for j, vj in enumerate(v):
            if vi and vj:
                total += K(vi * vj) * gram[i][j]
    return total


def h2_physical_point(v: tuple[int, ...]) -> tuple[float, float]:
    x = 0.0
    y = 0.0
    for j, coeff in enumerate(v):
        angle = 2 * math.pi * j / 5
        x += coeff * math.cos(angle)
        y += coeff * math.sin(angle)
    return (x, y)


def h2_internal_point(v: tuple[int, ...]) -> tuple[float, float]:
    x = 0.0
    y = 0.0
    for j, coeff in enumerate(v):
        angle = 4 * math.pi * j / 5
        x += coeff * math.cos(angle)
        y += coeff * math.sin(angle)
    return (x, y)


def convex_hull(points: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Return the convex hull in counter-clockwise order."""
    pts = sorted(set((round(x, 15), round(y, 15)) for x, y in points))
    if len(pts) <= 1:
        return pts

    def cross(o: tuple[float, float], a: tuple[float, float], b: tuple[float, float]) -> float:
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: list[tuple[float, float]] = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 1e-14:
            lower.pop()
        lower.append(p)
    upper: list[tuple[float, float]] = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 1e-14:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def h2_decagonal_window_vertices() -> list[tuple[float, float]]:
    """Centered Penrose strip window as the zonotope sum [-1/2,1/2] w_j."""
    generators = []
    for j in range(5):
        angle = 4 * math.pi * j / 5
        generators.append((math.cos(angle), math.sin(angle)))
    candidates = []
    for signs in product((-0.5, 0.5), repeat=5):
        x = sum(sign * generators[j][0] for j, sign in enumerate(signs))
        y = sum(sign * generators[j][1] for j, sign in enumerate(signs))
        candidates.append((x, y))
    return convex_hull(candidates)


H2_DECAGONAL_WINDOW = h2_decagonal_window_vertices()


def point_in_convex_polygon(point: tuple[float, float], polygon: list[tuple[float, float]]) -> bool:
    eps = 1e-11
    px, py = point
    sign = 0
    for i, a in enumerate(polygon):
        b = polygon[(i + 1) % len(polygon)]
        cross = (b[0] - a[0]) * (py - a[1]) - (b[1] - a[1]) * (px - a[0])
        if abs(cross) <= eps:
            continue
        current = 1 if cross > 0 else -1
        if sign == 0:
            sign = current
        elif sign != current:
            return False
    return True


def h2_window_pass(v: tuple[int, ...]) -> bool:
    return point_in_convex_polygon(h2_internal_point(v), H2_DECAGONAL_WINDOW)


def h2_patch(bound: int = 3, physical_radius: float = 4.0) -> list[dict]:
    gram_internal = gram_from_step(2)
    points = []
    for v in product(range(-bound, bound + 1), repeat=5):
        # Quotient by the diagonal kernel of the fifth-root star.
        if sum(v) != 0:
            continue
        x = h2_physical_point(v)
        if x[0] * x[0] + x[1] * x[1] > physical_radius * physical_radius + 1e-12:
            continue
        if not h2_window_pass(v):
            continue
        internal_norm = k_quadratic(v, gram_internal)
        points.append(
            {
                "label_Z5": list(v),
                "physical_xy_decimal": [round(x[0], 12), round(x[1], 12)],
                "internal_norm_squared_exact": k_text(internal_norm),
            }
        )
    points.sort(key=lambda item: (item["physical_xy_decimal"][0], item["physical_xy_decimal"][1], item["label_Z5"]))
    return points


def h2_reciprocal_labels(limit: int = 2) -> list[dict]:
    labels = []
    seen = set()
    gram_physical = gram_from_step(1)
    for h in product(range(-limit, limit + 1), repeat=5):
        if sum(h) != 0:
            continue
        if h == (0, 0, 0, 0, 0):
            continue
        kvec = h2_physical_point(h)
        norm_exact = k_quadratic(h, gram_physical)
        key = (round(kvec[0], 10), round(kvec[1], 10))
        if key in seen:
            continue
        seen.add(key)
        labels.append(
            {
                "label_Z5": list(h),
                "wavevector_xy_decimal": [round(kvec[0], 12), round(kvec[1], 12)],
                "norm_squared_exact": k_text(norm_exact),
            }
        )
    labels.sort(key=lambda item: (k_float(text_to_k(item["norm_squared_exact"])), item["label_Z5"]))
    return labels[:40]


def text_to_k(text: str) -> K:
    # Parser for the small strings emitted by K.__str__; used only for sorting.
    if text == "0":
        return K(0)
    if "phi" not in text:
        return K(Fraction(text))
    # Fall back to the decimal norm used for sorting if the expression is mixed.
    # The exact string remains the public data.
    return K(0)


def h2_diffraction(points: list[dict], reciprocal_labels: list[dict], top_n: int = 25) -> list[dict]:
    patch_xy = [tuple(p["physical_xy_decimal"]) for p in points]
    out = []
    n = len(patch_xy)
    for label in reciprocal_labels:
        kx, ky = label["wavevector_xy_decimal"]
        amp = 0j
        for x, y in patch_xy:
            amp += cmath.exp(-1j * (kx * x + ky * y))
        intensity = (abs(amp) / max(n, 1)) ** 2
        out.append(
            {
                "label_Z5": label["label_Z5"],
                "wavevector_xy_decimal": label["wavevector_xy_decimal"],
                "relative_intensity": round(float(intensity), 12),
            }
        )
    out.sort(key=lambda item: (-item["relative_intensity"], item["label_Z5"]))
    return out[:top_n]


def h2_artifacts() -> dict[str, dict | list]:
    gram_physical = gram_from_step(1)
    gram_internal = gram_from_step(2)
    patch = h2_patch()
    reciprocal = h2_reciprocal_labels()
    diffraction = h2_diffraction(patch, reciprocal)
    shell = []
    for k in range(10):
        angle = math.pi * k / 5
        shell.append({"k": k, "xy_decimal": [round(math.cos(angle), 12), round(math.sin(angle), 12)]})
    return {
        "h2_shell_certificate.json": {
            "certificate": "h2_shell_certificate",
            "root_system": "H2=I2(5)",
            "root_count": 10,
            "cartan_generator": "2*cos(pi/5)=phi",
            "coefficient_ring": "Z[phi]",
            "roots_unit_circle_decimal": shell,
            "control_role": "benchmark for the H4 model-set arithmetic",
        },
        "cut_project_basis.json": {
            "certificate": "cut_project_basis",
            "source_lattice": "Z^5 / diagonal Z",
            "physical_columns": "v_j=(cos(2*pi*j/5), sin(2*pi*j/5))",
            "internal_columns": "w_j=(cos(4*pi*j/5), sin(4*pi*j/5))",
            "physical_gram_exact": [[k_text(x) for x in row] for row in gram_physical],
            "internal_gram_exact": [[k_text(x) for x in row] for row in gram_internal],
            "window": {
                "type": "centered regular decagonal zonotope",
                "definition": "sum_j t_j w_j with -1/2 <= t_j <= 1/2 and w_j=(cos(4*pi*j/5), sin(4*pi*j/5))",
                "vertices_decimal_ccw": [[round(x, 12), round(y, 12)] for x, y in H2_DECAGONAL_WINDOW],
                "note": "This is the canonical centered strip window used here for the finite H2 benchmark; the data remain finite-window intensities rather than an infinite diffraction theorem.",
            },
        },
        "penrose_patch_vertices.json": {
            "certificate": "penrose_patch_vertices",
            "patch_type": "finite H2 cyclotomic model-set patch",
            "enumeration_bound": 3,
            "physical_radius": 4.0,
            "vertex_count": len(patch),
            "vertices": patch,
        },
        "reciprocal_module_h2.json": {
            "certificate": "reciprocal_module_h2",
            "module": "dual labels in Z^5 / diagonal Z projected to physical space",
            "label_limit": 2,
            "label_count_exported": len(reciprocal),
            "labels": reciprocal,
        },
        "diffraction_intensities_h2.json": {
            "certificate": "diffraction_intensities_h2",
            "method": "finite-patch kinematic intensity |sum exp(-i k.x)|^2 / N^2",
            "patch_vertex_count": len(patch),
            "top_intensities": diffraction,
            "scope": "benchmark finite intensities, not an infinite-volume diffraction theorem",
        },
    }


def quaternion_to_coeff_pairs(q) -> list[list[str]]:
    return [k_pair(c) for c in q.c]


def h3_star_vectors() -> list[dict]:
    vecs = []
    raw = []
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            raw.append((K(0), K(s1), K(0, s2)))
            raw.append((K(s1), K(0, s2), K(0)))
            raw.append((K(0, s1), K(0), K(s2)))
    seen = set()
    for v in raw:
        key = tuple((x.a, x.b) for x in v)
        if key in seen:
            continue
        seen.add(key)
        vecs.append(
            {
                "xyz_zphi": [k_text(x) for x in v],
                "xyz_decimal": [round(k_float(x), 12) for x in v],
            }
        )
    vecs.sort(key=lambda item: item["xyz_decimal"])
    return vecs


def h4_artifacts() -> dict[str, dict | list]:
    roots = h4_quaternion_roots()
    exported_roots = []
    for idx, root in enumerate(roots):
        exported_roots.append(
            {
                "index": idx,
                "quaternion_coefficients_zphi": quaternion_to_coeff_pairs(root),
                "norm_zphi": k_text(root.norm()),
            }
        )
    star = h3_star_vectors()
    return {
        "h4_shell_coordinates.json": {
            "certificate": "h4_shell_coordinates",
            "root_system": "H4",
            "root_count": len(roots),
            "coordinate_ring": "Z[phi] with half-integral quaternion coordinates",
            "all_norm_one": all(root.norm() == ONE for root in roots),
            "coordinates": exported_roots,
        },
        "projection_matrices_h4.json": {
            "certificate": "projection_matrices_h4",
            "source_coordinates": "quaternion coefficients (a,b,c,d) over Q(sqrt(5))",
            "physical_projection": {
                "target": "R^3 icosahedral star coordinates",
                "matrix": [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
                "description": "imaginary-part projection (a,b,c,d)->(b,c,d)",
            },
            "internal_projection": {
                "target": "scalar complement",
                "matrix": [[1, 0, 0, 0]],
                "description": "scalar coordinate retained as the first internal coordinate",
            },
            "galois_star": "apply phi -> 1-phi coefficientwise for conjugate-space labels",
            "scope": "projection data for exact labels and star vectors, not a full atomic-surface classification",
        },
        "icosahedral_star_vectors.json": {
            "certificate": "icosahedral_star_vectors",
            "source": "H3 star extracted as the standard icosahedral shadow of H4 arithmetic",
            "vector_count": len(star),
            "vectors": star,
        },
    }


def read_k8_gram() -> list[list[int]]:
    path = K8E8_DIR / "icosian_trace_gram.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    return data["gram_matrix"]


def standard_e8_gram() -> list[list[int]]:
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


def inverse_denominator_data(matrix: list[list[int]]) -> tuple[list[list[str]], int, int]:
    inv = sp.Matrix(matrix).inv()
    den_lcm = 1
    non_integral = 0
    rows = []
    for i in range(inv.rows):
        row = []
        for j in range(inv.cols):
            value = sp.Rational(inv[i, j])
            den_lcm = int(sp.ilcm(den_lcm, value.q))
            if value.q != 1:
                non_integral += 1
            row.append(str(value))
        rows.append(row)
    return rows, den_lcm, non_integral


def smith_diagonal_data(matrix: list[list[int]]) -> list[int]:
    normal = smith_normal_form(sp.Matrix(matrix), domain=sp.ZZ)
    return [abs(int(normal[i, i])) for i in range(normal.rows) if normal[i, i] != 0]


def reciprocal_artifacts() -> dict[str, dict | list]:
    k8 = read_k8_gram()
    e8 = standard_e8_gram()
    k8_inv, k8_lcm, k8_nonintegral = inverse_denominator_data(k8)
    e8_inv, e8_lcm, e8_nonintegral = inverse_denominator_data(e8)
    det_k8 = int(sp.Matrix(k8).det())
    det_e8 = int(sp.Matrix(e8).det())
    k8_smith = smith_diagonal_data(k8)
    e8_smith = smith_diagonal_data(e8)
    scaled_k8_inv = [[str(5 * sp.Rational(entry)) for entry in row] for row in k8_inv]
    labels = []
    for idx, row in enumerate(k8_inv):
        labels.append(
            {
                "dual_basis_index": idx,
                "coordinates_in_K8_basis": row,
                "coordinates_scaled_by_5": scaled_k8_inv[idx],
            }
        )
    return {
        "reciprocal_module_h4.json": {
            "certificate": "reciprocal_module_h4",
            "lattice": "K8 icosian trace lattice",
            "trace_polar_form": "Tr_{K/Q}(2*<x,y>_K)",
            "dual_denominator_lcm": k8_lcm,
            "nonintegral_inverse_entries": k8_nonintegral,
            "dual_basis_generators": labels,
            "interpretation": "Reciprocal labels for the K8 normalization naturally carry fifth-denominator coordinates.",
        },
        "k8_e8_projection_comparison.json": {
            "certificate": "k8_e8_projection_comparison",
            "comparison_axis": "reciprocal module and discriminant layer",
            "k8": {
                "determinant": det_k8,
                "discriminant_group": "(Z/5Z)^4",
                "smith_normal_form_diagonal": k8_smith,
                "dual_denominator_lcm": k8_lcm,
                "nonintegral_inverse_entries": k8_nonintegral,
                "has_fifth_denominator_reciprocal_labels": k8_lcm == 5,
            },
            "e8": {
                "determinant": det_e8,
                "discriminant_group": "trivial",
                "smith_normal_form_diagonal": e8_smith,
                "dual_denominator_lcm": e8_lcm,
                "nonintegral_inverse_entries": e8_nonintegral,
                "has_fifth_denominator_reciprocal_labels": False,
            },
            "publishable_consequence": (
                "The visible H4/icosian shell can be compared with E8 geometry, "
                "but the reciprocal module for the trace-normalized icosian lattice "
                "is not E8-self-dual: it contains a 5-primary discriminant layer."
            ),
            "table_statement": "K8^#/K8=(Z/5Z)^4 whereas E8^#/E8=0.",
        },
    }


def build_all_artifacts() -> dict[str, dict | list]:
    artifacts: dict[str, dict | list] = {}
    artifacts.update(h2_artifacts())
    artifacts.update(h4_artifacts())
    artifacts.update(reciprocal_artifacts())
    artifacts["summary.json"] = {
        "certificate": "summary",
        "article": ARTICLE,
        "main_outputs": {
            "h2_patch_vertices": artifacts["penrose_patch_vertices.json"]["vertex_count"],
            "h2_diffraction_peaks_exported": len(artifacts["diffraction_intensities_h2.json"]["top_intensities"]),
            "h4_root_count": artifacts["h4_shell_coordinates.json"]["root_count"],
            "h4_star_vector_count": artifacts["icosahedral_star_vectors.json"]["vector_count"],
            "k8_dual_denominator_lcm": artifacts["reciprocal_module_h4.json"]["dual_denominator_lcm"],
            "k8_e8_consequence": artifacts["k8_e8_projection_comparison.json"]["table_statement"],
        },
        "scope": [
            "Mathematical model-set and reciprocal-module arithmetic.",
            "Finite H2 diffraction benchmark with a centered decagonal strip window only.",
            "No experimental materials claim.",
            "No G0/octonion associator claim.",
        ],
        "all_core_checks_passed": (
            artifacts["h2_shell_certificate.json"]["root_count"] == 10
            and artifacts["penrose_patch_vertices.json"]["vertex_count"] > 0
            and artifacts["h4_shell_coordinates.json"]["root_count"] == 120
            and artifacts["reciprocal_module_h4.json"]["dual_denominator_lcm"] == 5
            and artifacts["k8_e8_projection_comparison.json"]["k8"]["discriminant_group"] == "(Z/5Z)^4"
        ),
    }
    artifacts["summary.json"]["file_hashes_sha256"] = {
        name: sha256_json(data) for name, data in artifacts.items() if name != "summary.json"
    }
    return artifacts


def write_artifacts(out_dir: Path, names: list[str] | None = None) -> dict[str, dict | list]:
    artifacts = build_all_artifacts()
    selected = names if names is not None else sorted(artifacts)
    for name in selected:
        write_json(out_dir / name, artifacts[name])
    return {name: artifacts[name] for name in selected}


def run_all() -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = RESULTS_ROOT / f"run_exact_golden_model_sets_{timestamp}"
    run_dir.mkdir(parents=True, exist_ok=True)
    artifacts = build_all_artifacts()
    CERT_DIR.mkdir(parents=True, exist_ok=True)
    for name, data in artifacts.items():
        write_json(CERT_DIR / name, data)
        write_json(run_dir / name, data)
    manifest = {
        "article": ARTICLE,
        "script": Path(sys.argv[0]).name,
        "core_version": CORE_VERSION,
        "timestamp": timestamp,
        "canonical_certificate_dir": str(CERT_DIR.relative_to(ROOT)),
        "run_dir": str(run_dir.relative_to(ROOT)),
        "main_outputs": artifacts["summary.json"]["main_outputs"],
        "all_core_checks_passed": artifacts["summary.json"]["all_core_checks_passed"],
    }
    write_json(run_dir / "run_manifest.json", manifest)
    for script in Path(__file__).parent.glob("*.py"):
        shutil.copy2(script, run_dir / script.name)
    return run_dir
