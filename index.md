---
title: Non-crystallographic integers — Certificates
---

# Non-crystallographic integers — Certificate archive

Exact-arithmetic verification scripts and machine-readable certificates
accompanying the article

> **D. Corradetti**, *Non-crystallographic systems of integers over composition
> algebras*, 2026.

This site is the public supplementary archive referenced in the *Computational
methods* section of the paper. All numerical claims in the article — finite
shell cardinalities, Gram matrix determinants, projective-line tallies,
denominator searches, trace-polar discriminant quotients — are reproduced here
in exact arithmetic over `Z[φ]`, with one JSON certificate per logical step.

## Contents

- [Verification protocol (P1)–(P6)](docs/protocol.md)
- [Paper ↔ certificate correspondence](docs/correspondence.md)
- [Source code](scripts/)
- [JSON certificates](certificates/)
- [Companion supplement: exact golden model sets (Article 17)](supplement/exact_golden_model_sets/)
- [Source repository on GitHub](https://github.com/DCorradetti/Non-crystallographic-Integers)

## Quick verification

Clone the repository and run, from its root:

```bash
python verify.py
```

The implementation has no third-party dependencies; only the Python standard
library (≥ 3.11) is required. Total runtime: under one minute on a standard
laptop.

A SHA-256 digest of the certificate tree is printed at the end and stored in
[`certificates/SHA256SUMS`](certificates/SHA256SUMS).

## Headline results

- **Icosian double `G_0 = I ⊕ Iℓ`** is a free `Z[φ]`-order of rank 8 in
  `O(K)`, with a 240-element finite shell of type `H_4 ⊕ H_4`, closed under
  reflections, with Cartan coefficients in `Z[φ]`, and genuinely octonionic
  (nonzero associator `[i,j,ℓ] = 2kℓ`).
  → [`certificates/icosian_double/`](certificates/icosian_double/)
- **Self-duality `G_0^# = G_0`** with respect to the polar norm pairing;
  `det G_I = φ²`, `det G_{G_0} = φ⁴ = 2 + 3φ`, of field norm `1`.
  → [`certificates/g3/dual_discriminant_certificate.json`](certificates/g3/dual_discriminant_certificate.json)
- **No denominator-2 G3 candidate**: all `21845` projective lines of
  `(Z[φ]/2)^8 ≅ F_4^8` fail one of three structural filters
  (170 not-mixed, 16320 conjugation-stable failure, 5355 integrally unpaired).
  → [`certificates/g3/denominator2_all_subspaces_no_go.json`](certificates/g3/denominator2_all_subspaces_no_go.json)
- **No ramified `√5`-denominator G3 candidate**: all `97656` projective lines
  of `(Z[φ]/√5)^8 ≅ F_5^8` fail the polar-pairing filter; the structural
  reason is non-degeneracy of `B mod √5`.
  → [`certificates/g3/denominator_sqrt5_all_subspaces_no_go.json`](certificates/g3/denominator_sqrt5_all_subspaces_no_go.json)
- **No G3-B' candidate inside the icosian-double discriminant tower**: the
  trace-polar quotient `G_{0,Z}^# / G_0 ≅ (Z/5Z)^8` carries a split form of
  type `O^+(8,5)` with `19656` isotropic projective lines, every one of which
  generates a full-dimension stable closure.
  → [`certificates/g3_b_prime/`](certificates/g3_b_prime/)

## Companion supplements

The [`supplement/`](supplement/) tree hosts self-contained reproducibility
packages for companion articles.

- **[`exact_golden_model_sets/`](supplement/exact_golden_model_sets/)** —
  D. Corradetti, *Exact golden arithmetic of H2 and H4 model sets*. Bundles
  v2 scripts, JSON certificates, and Markdown audit notes. Requires `sympy`
  in addition to the Python standard library.

## License

- Source code: MIT
- Certificates and documentation: CC-BY-4.0
