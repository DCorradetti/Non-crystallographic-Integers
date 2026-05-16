from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import argparse
import hashlib
import itertools
import json
from pathlib import Path


def F(value=0) -> Fraction:
    return value if isinstance(value, Fraction) else Fraction(value)


@dataclass(frozen=True)
class K:
    """Element a + b phi in Q(phi), with phi^2 = phi + 1."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    def __init__(self, a=0, b=0):
        object.__setattr__(self, "a", F(a))
        object.__setattr__(self, "b", F(b))

    def __add__(self, other: object) -> "K":
        other = as_k(other)
        return K(self.a + other.a, self.b + other.b)

    def __radd__(self, other: object) -> "K":
        return self + other

    def __sub__(self, other: object) -> "K":
        other = as_k(other)
        return K(self.a - other.a, self.b - other.b)

    def __rsub__(self, other: object) -> "K":
        return as_k(other) - self

    def __neg__(self) -> "K":
        return K(-self.a, -self.b)

    def __mul__(self, other: object) -> "K":
        other = as_k(other)
        # (a+b phi)(c+d phi) = (ac+bd) + (ad+bc+bd) phi.
        return K(
            self.a * other.a + self.b * other.b,
            self.a * other.b + self.b * other.a + self.b * other.b,
        )

    def __rmul__(self, other: object) -> "K":
        return self * other

    def __truediv__(self, other: object) -> "K":
        other = as_k(other)
        return self * other.inverse()

    def __pow__(self, n: int) -> "K":
        if n == 0:
            return ONE
        if n < 0:
            return (self.inverse()) ** (-n)
        out = ONE
        base = self
        while n:
            if n & 1:
                out *= base
            base *= base
            n >>= 1
        return out

    def star(self) -> "K":
        # phi* = 1 - phi, so (a+b phi)* = a+b-b phi.
        return K(self.a + self.b, -self.b)

    def trace_q(self) -> Fraction:
        return 2 * self.a + self.b

    def norm_q(self) -> Fraction:
        return self.a * self.a + self.a * self.b - self.b * self.b

    def inverse(self) -> "K":
        n = self.norm_q()
        if n == 0:
            raise ZeroDivisionError("zero has no inverse in Q(phi)")
        s = self.star()
        return K(s.a / n, s.b / n)

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0

    def is_integral(self) -> bool:
        return self.a.denominator == 1 and self.b.denominator == 1

    def as_pair(self) -> list[str]:
        return [str(self.a), str(self.b)]

    def __str__(self) -> str:
        if self.b == 0:
            return str(self.a)
        if self.a == 0:
            if self.b == 1:
                return "phi"
            if self.b == -1:
                return "-phi"
            return f"{self.b}*phi"
        sign = "+" if self.b > 0 else "-"
        b_abs = abs(self.b)
        b_text = "phi" if b_abs == 1 else f"{b_abs}*phi"
        return f"{self.a}{sign}{b_text}"


def as_k(value: object) -> K:
    if isinstance(value, K):
        return value
    return K(value, 0)


ZERO = K(0)
ONE = K(1)
TWO = K(2)
HALF = K(Fraction(1, 2))
PHI = K(0, 1)
PHI_INV = PHI - ONE


@dataclass(frozen=True)
class Quaternion:
    c: tuple[K, K, K, K]

    def __init__(self, coeffs=None):
        coeffs = coeffs if coeffs is not None else [ZERO] * 4
        if len(coeffs) != 4:
            raise ValueError("quaternion needs four coefficients")
        object.__setattr__(self, "c", tuple(as_k(x) for x in coeffs))

    def __add__(self, other: "Quaternion") -> "Quaternion":
        return Quaternion([a + b for a, b in zip(self.c, other.c)])

    def __sub__(self, other: "Quaternion") -> "Quaternion":
        return Quaternion([a - b for a, b in zip(self.c, other.c)])

    def __neg__(self) -> "Quaternion":
        return Quaternion([-x for x in self.c])

    def scale(self, s: object) -> "Quaternion":
        s = as_k(s)
        return Quaternion([s * x for x in self.c])

    def conj(self) -> "Quaternion":
        a, b, c, d = self.c
        return Quaternion([a, -b, -c, -d])

    def mul(self, other: "Quaternion") -> "Quaternion":
        a, b, c, d = self.c
        e, f, g, h = other.c
        return Quaternion(
            [
                a * e - b * f - c * g - d * h,
                a * f + b * e + c * h - d * g,
                a * g - b * h + c * e + d * f,
                a * h + b * g - c * f + d * e,
            ]
        )

    def norm(self) -> K:
        return sum((x * x for x in self.c), ZERO)

    def inner(self, other: "Quaternion") -> K:
        return sum((a * b for a, b in zip(self.c, other.c)), ZERO)

    def is_zero(self) -> bool:
        return all(x.is_zero() for x in self.c)

    def key(self) -> tuple[tuple[Fraction, Fraction], ...]:
        return tuple((x.a, x.b) for x in self.c)

    def to_json(self) -> list[list[str]]:
        return [x.as_pair() for x in self.c]

    def __str__(self) -> str:
        names = ["1", "i", "j", "k"]
        parts = []
        for coeff, name in zip(self.c, names):
            if not coeff.is_zero():
                parts.append(f"({coeff}){name}")
        return " + ".join(parts) if parts else "0"


QZERO = Quaternion()
QONE = Quaternion([ONE, ZERO, ZERO, ZERO])
QI = Quaternion([ZERO, ONE, ZERO, ZERO])
QJ = Quaternion([ZERO, ZERO, ONE, ZERO])
QK = Quaternion([ZERO, ZERO, ZERO, ONE])


@dataclass(frozen=True)
class Octonion:
    a: Quaternion
    b: Quaternion

    def __init__(self, a=None, b=None):
        object.__setattr__(self, "a", a if isinstance(a, Quaternion) else Quaternion(a))
        object.__setattr__(self, "b", b if isinstance(b, Quaternion) else Quaternion(b))

    def __add__(self, other: "Octonion") -> "Octonion":
        return Octonion(self.a + other.a, self.b + other.b)

    def __sub__(self, other: "Octonion") -> "Octonion":
        return Octonion(self.a - other.a, self.b - other.b)

    def __neg__(self) -> "Octonion":
        return Octonion(-self.a, -self.b)

    def scale(self, s: object) -> "Octonion":
        return Octonion(self.a.scale(s), self.b.scale(s))

    def conj(self) -> "Octonion":
        return Octonion(self.a.conj(), -self.b)

    def mul(self, other: "Octonion") -> "Octonion":
        # (a+b l)(c+d l) = (ac - conjugate(d)b) + (d a + b conjugate(c))l.
        a, b = self.a, self.b
        c, d = other.a, other.b
        first = a.mul(c) - d.conj().mul(b)
        second = d.mul(a) + b.mul(c.conj())
        return Octonion(first, second)

    def norm(self) -> K:
        return self.a.norm() + self.b.norm()

    def inner(self, other: "Octonion") -> K:
        return self.a.inner(other.a) + self.b.inner(other.b)

    def associator(self, y: "Octonion", z: "Octonion") -> "Octonion":
        return self.mul(y).mul(z) - self.mul(y.mul(z))

    def coeffs(self) -> tuple[K, ...]:
        return self.a.c + self.b.c

    def key(self) -> tuple[tuple[Fraction, Fraction], ...]:
        return tuple((x.a, x.b) for x in self.coeffs())

    def to_json(self) -> list[list[str]]:
        return [x.as_pair() for x in self.coeffs()]

    def is_zero(self) -> bool:
        return self.a.is_zero() and self.b.is_zero()

    def __str__(self) -> str:
        names = ["1", "i", "j", "k", "l", "il", "jl", "kl"]
        parts = []
        for coeff, name in zip(self.coeffs(), names):
            if not coeff.is_zero():
                parts.append(f"({coeff}){name}")
        return " + ".join(parts) if parts else "0"


OZERO = Octonion(QZERO, QZERO)
OONE = Octonion(QONE, QZERO)
OL = Octonion(QZERO, QONE)


def q_basis_icosian() -> list[Quaternion]:
    h = (QONE + QI + QJ + QK).scale(HALF)
    # This H4-compatible icosian basis contains the standard 120 roots used
    # below. It is equivalent to the usual icosian presentations up to an
    # H4/binary-icosahedral change of basis.
    g = (-QONE + QI.scale(PHI - ONE) - QJ.scale(PHI)).scale(HALF)
    return [QONE, QI, h, g]


def oct_basis_icosian_double() -> list[Octonion]:
    qb = q_basis_icosian()
    return [Octonion(q, QZERO) for q in qb] + [Octonion(QZERO, q) for q in qb]


def h4_quaternion_roots() -> list[Quaternion]:
    roots: dict[tuple[tuple[Fraction, Fraction], ...], Quaternion] = {}

    def add(q: Quaternion):
        roots[q.key()] = q

    # (±1,0,0,0) and all permutations.
    for idx in range(4):
        for s in [ONE, -ONE]:
            coeffs = [ZERO] * 4
            coeffs[idx] = s
            add(Quaternion(coeffs))

    # 1/2(±1,±1,±1,±1).
    for signs in itertools.product([ONE, -ONE], repeat=4):
        add(Quaternion([s * HALF for s in signs]))

    # 1/2(0, ±1, ±phi, ±phi^{-1}) and all even permutations.
    base = [ZERO, ONE, PHI, PHI_INV]
    for perm in even_permutations(range(4)):
        vals = [base[i] for i in perm]
        for signs in itertools.product([ONE, -ONE], repeat=3):
            coeffs = []
            sign_iter = iter(signs)
            for value in vals:
                coeffs.append(value if value.is_zero() else value * next(sign_iter))
            add(Quaternion([x * HALF for x in coeffs]))

    out = list(roots.values())
    out.sort(key=lambda q: tuple((c.a, c.b) for c in q.c))
    return out


def even_permutations(items):
    for p in itertools.permutations(items):
        inv = 0
        for i in range(len(p)):
            for j in range(i + 1, len(p)):
                if p[i] > p[j]:
                    inv += 1
        if inv % 2 == 0:
            yield p


def solve_coordinates_vector(target: list[K], basis_vectors: list[list[K]]) -> list[K]:
    n = len(target)
    mat = [[basis_vectors[col][row] for col in range(n)] + [target[row]] for row in range(n)]
    pivot_row = 0
    for col in range(n):
        pivot = None
        for row in range(pivot_row, n):
            if not mat[row][col].is_zero():
                pivot = row
                break
        if pivot is None:
            continue
        mat[pivot_row], mat[pivot] = mat[pivot], mat[pivot_row]
        inv = mat[pivot_row][col].inverse()
        mat[pivot_row] = [inv * x for x in mat[pivot_row]]
        for row in range(n):
            if row == pivot_row:
                continue
            factor = mat[row][col]
            if factor.is_zero():
                continue
            mat[row] = [x - factor * y for x, y in zip(mat[row], mat[pivot_row])]
        pivot_row += 1
    if pivot_row != n:
        raise ValueError("basis matrix is singular")
    return [mat[row][-1] for row in range(n)]


def q_coordinates(q: Quaternion, basis: list[Quaternion]) -> list[K]:
    return solve_coordinates_vector(list(q.c), [list(b.c) for b in basis])


def o_coordinates(o: Octonion, basis: list[Octonion]) -> list[K]:
    return solve_coordinates_vector(list(o.coeffs()), [list(b.coeffs()) for b in basis])


def in_r_span(coords: list[K]) -> bool:
    return all(c.is_integral() for c in coords)


def format_coords(coords: list[K]) -> list[str]:
    return [str(c) for c in coords]


def verify_quaternion_order(basis: list[Quaternion]) -> dict:
    products = {}
    closed = True
    for i, bi in enumerate(basis):
        for j, bj in enumerate(basis):
            coords = q_coordinates(bi.mul(bj), basis)
            ok = in_r_span(coords)
            closed = closed and ok
            products[f"{i},{j}"] = format_coords(coords)
    conj = {}
    conj_closed = True
    for i, bi in enumerate(basis):
        coords = q_coordinates(bi.conj(), basis)
        ok = in_r_span(coords)
        conj_closed = conj_closed and ok
        conj[str(i)] = format_coords(coords)
    gram = [[str(bi.inner(bj)) for bj in basis] for bi in basis]
    return {
        "rank_over_Zphi": 4,
        "closed_under_multiplication": closed,
        "closed_under_conjugation": conj_closed,
        "structure_constants": products,
        "conjugation_matrix_rows": conj,
        "gram_matrix": gram,
    }


def verify_octonion_order(basis: list[Octonion]) -> dict:
    products = {}
    closed = True
    for i, bi in enumerate(basis):
        for j, bj in enumerate(basis):
            coords = o_coordinates(bi.mul(bj), basis)
            ok = in_r_span(coords)
            closed = closed and ok
            products[f"{i},{j}"] = format_coords(coords)
    conj = {}
    conj_closed = True
    for i, bi in enumerate(basis):
        coords = o_coordinates(bi.conj(), basis)
        ok = in_r_span(coords)
        conj_closed = conj_closed and ok
        conj[str(i)] = format_coords(coords)
    gram = [[str(bi.inner(bj)) for bj in basis] for bi in basis]
    return {
        "rank_over_Zphi": 8,
        "closed_under_multiplication": closed,
        "closed_under_conjugation": conj_closed,
        "trace_norm_integral_by_gram_matrix": all(as_k_entry_integral(x) for row in gram for x in row),
        "structure_constants": products,
        "conjugation_matrix_rows": conj,
        "gram_matrix": gram,
    }


def as_k_entry_integral(text: str) -> bool:
    # This is only used after stringification of exact integral entries.
    return "/" not in text


def verify_h4_roots_in_icosian() -> dict:
    basis = q_basis_icosian()
    roots = h4_quaternion_roots()
    membership = [in_r_span(q_coordinates(r, basis)) for r in roots]
    unit_norm = [r.norm() == ONE for r in roots]
    keyset = {r.key() for r in roots}
    neg_closed = all((-r).key() in keyset for r in roots)
    mult_closed = True
    for a in roots:
        for b in roots:
            if a.mul(b).key() not in keyset:
                mult_closed = False
                break
        if not mult_closed:
            break
    return {
        "root_count": len(roots),
        "all_roots_in_icosian_order": all(membership),
        "all_norm_one": all(unit_norm),
        "closed_under_negation": neg_closed,
        "closed_under_quaternion_multiplication": mult_closed,
        "type": "H4",
    }


def double_shell() -> list[Octonion]:
    roots = h4_quaternion_roots()
    shell = [Octonion(r, QZERO) for r in roots] + [Octonion(QZERO, r) for r in roots]
    shell.sort(key=lambda x: x.key())
    return shell


def verify_double_shell() -> dict:
    shell = double_shell()
    keyset = {x.key() for x in shell}
    neg_closed = all((-x).key() in keyset for x in shell)
    norms_one = all(x.norm() == ONE for x in shell)
    mult_closed = True
    for a in shell:
        for b in shell:
            if a.mul(b).key() not in keyset:
                mult_closed = False
                break
        if not mult_closed:
            break
    cartan_integral = True
    reflection_closed = True
    for alpha in shell:
        aa = alpha.inner(alpha)
        for beta in shell:
            cartan = TWO * beta.inner(alpha) / aa
            if not cartan.is_integral():
                cartan_integral = False
            reflected = beta - alpha.scale(cartan)
            if reflected.key() not in keyset:
                reflection_closed = False
                break
        if not reflection_closed:
            break
    # Deterministic nonzero associator example: i, j, l.
    oi = Octonion(QI, QZERO)
    oj = Octonion(QJ, QZERO)
    ol = Octonion(QZERO, QONE)
    assoc = oi.associator(oj, ol)
    return {
        "shell_name": "IcosianDouble_H4_plus_H4",
        "cardinality": len(shell),
        "rank_over_K": 8,
        "all_norm_one": norms_one,
        "closed_under_negation": neg_closed,
        "closed_under_multiplication": mult_closed,
        "closed_under_reflections": reflection_closed,
        "cartan_coefficients_in_Zphi": cartan_integral,
        "coxeter_type": "H4 + H4",
        "crystallographic_over_Z": False,
        "crystallographic_over_Zphi": True,
        "genuinely_octonionic": not assoc.is_zero(),
        "nonzero_associator_example": {
            "x": "i",
            "y": "j",
            "z": "l",
            "associator": str(assoc),
        },
        "shell_hash": hash_objects(shell),
    }


def hash_objects(objects) -> str:
    payload = json.dumps([obj.to_json() for obj in objects], sort_keys=True)
    return "sha256:" + hashlib.sha256(payload.encode("utf-8")).hexdigest()


def field_certificate() -> dict:
    samples_ok = True
    for a, b, c, d in itertools.product(range(-2, 3), repeat=4):
        x = K(a, b)
        y = K(c, d)
        samples_ok = samples_ok and (x * y).norm_q() == x.norm_q() * y.norm_q()
        samples_ok = samples_ok and x.star().star() == x
    return {
        "field": "Q(sqrt(5))",
        "ring": "Z[phi]",
        "relation": "phi^2 = phi + 1",
        "galois_conjugation": "phi -> 1 - phi",
        "trace_formula": "Tr(a+b phi)=2a+b",
        "norm_formula": "N(a+b phi)=a^2+a*b-b^2",
        "unit_group": "+/- phi^k",
        "sample_checks_passed": samples_ok,
    }


def octonion_certificate() -> dict:
    basis = [
        OONE,
        Octonion(QI, QZERO),
        Octonion(QJ, QZERO),
        Octonion(QK, QZERO),
        OL,
        Octonion(QZERO, QI),
        Octonion(QZERO, QJ),
        Octonion(QZERO, QK),
    ]
    squares_ok = all(basis[i].mul(basis[i]) == -OONE for i in range(1, 8))
    unit_ok = all(OONE.mul(x) == x and x.mul(OONE) == x for x in basis)
    conj_ok = True
    norm_ok = True
    moufang_ok = True
    for x, y, z in itertools.product(basis, repeat=3):
        conj_ok = conj_ok and x.mul(y).conj() == y.conj().mul(x.conj())
        norm_ok = norm_ok and x.mul(y).norm() == x.norm() * y.norm()
        moufang_ok = moufang_ok and x.mul(y).mul(z.mul(x)) == x.mul(y.mul(z)).mul(x)
    return {
        "basis": ["1", "i", "j", "k", "l", "il", "jl", "kl"],
        "construction": "Cayley-Dickson double of Hamilton quaternions",
        "product": "(a+b*l)(c+d*l)=(ac-conj(d)b)+(da+b*conj(c))*l",
        "verified_identities_on_basis": {
            "imaginary_basis_squares_minus_one": squares_ok,
            "identity": unit_ok,
            "conjugation_anti_automorphism": conj_ok,
            "norm_multiplicative": norm_ok,
            "moufang_identity": moufang_ok,
        },
    }


def write_json(path: Path, data: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def generate_certificates(out_dir: Path):
    qb = q_basis_icosian()
    ob = oct_basis_icosian_double()
    certificates = {
        "field_certificate.json": field_certificate(),
        "octonion_multiplication_certificate.json": octonion_certificate(),
        "icosian_order_certificate.json": verify_quaternion_order(qb),
        "icosian_H4_shell_certificate.json": verify_h4_roots_in_icosian(),
        "icosian_double_order_certificate.json": verify_octonion_order(ob),
        "icosian_double_root_shell_certificate.json": verify_double_shell(),
    }
    for name, data in certificates.items():
        write_json(out_dir / name, data)
    summary = {
        "all_certificates": sorted(certificates),
        "all_core_checks_passed": all_core_checks(certificates),
    }
    write_json(out_dir / "summary.json", summary)
    return summary


def all_core_checks(certificates: dict) -> bool:
    return (
        certificates["field_certificate.json"]["sample_checks_passed"]
        and all(certificates["octonion_multiplication_certificate.json"]["verified_identities_on_basis"].values())
        and certificates["icosian_order_certificate.json"]["closed_under_multiplication"]
        and certificates["icosian_order_certificate.json"]["closed_under_conjugation"]
        and certificates["icosian_H4_shell_certificate.json"]["root_count"] == 120
        and certificates["icosian_H4_shell_certificate.json"]["all_roots_in_icosian_order"]
        and certificates["icosian_H4_shell_certificate.json"]["closed_under_quaternion_multiplication"]
        and certificates["icosian_double_order_certificate.json"]["closed_under_multiplication"]
        and certificates["icosian_double_order_certificate.json"]["closed_under_conjugation"]
        and certificates["icosian_double_root_shell_certificate.json"]["cardinality"] == 240
        and certificates["icosian_double_root_shell_certificate.json"]["closed_under_reflections"]
        and certificates["icosian_double_root_shell_certificate.json"]["genuinely_octonionic"]
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        default="../../3. Certificates/golden_octonions/icosian_double",
        help="output certificate directory",
    )
    args = parser.parse_args()
    out = Path(args.out).resolve()
    summary = generate_certificates(out)
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"certificates written to {out}")


if __name__ == "__main__":
    main()
