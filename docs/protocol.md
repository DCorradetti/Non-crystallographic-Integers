---
title: Verification protocol (P1)–(P6)
---

# Verification protocol (P1)–(P6)

This page reproduces the protocol from §8 *Computational methods* of the
article. Each step lists its mathematical content, the script(s) that
execute it, and the JSON certificate(s) it produces.

All arithmetic is exact:
- elements of `Z[φ]` are represented as integer pairs `(a, b)` standing for
  `a + b φ`, with multiplication reduced via `φ² = φ + 1`;
- elements of `K = Q(√5)` are represented as pairs of `fractions.Fraction`;
- the icosian ring `I` uses the four-element basis
  `(1, i, (1+i+j+k)/2, (-1+(φ-1)i-φj)/2)`;
- the icosian double `G_0 = I ⊕ Iℓ` uses the 8-element basis obtained by
  appending `ℓ`-multiples of the above; multiplication is via the
  Cayley–Dickson formula `(a + bℓ)(c + dℓ) = (ac - d̄b) + (da + b c̄)ℓ`.

## (P1) Closure and integrality

**Content.** `I` and `G_0` are closed under multiplication and conjugation;
`Tr_A` and `Norm_A` are integral on every basis pair.

**Script.** [`scripts/golden_octonions.py`](../scripts/golden_octonions.py)

**Certificates.** [`certificates/icosian_double/icosian_order_certificate.json`](../certificates/icosian_double/icosian_order_certificate.json),
[`certificates/icosian_double/icosian_double_order_certificate.json`](../certificates/icosian_double/icosian_double_order_certificate.json),
[`certificates/icosian_double/octonion_multiplication_certificate.json`](../certificates/icosian_double/octonion_multiplication_certificate.json)

## (P2) Enumeration of the `H_4` shell

**Content.** The 120-element shell `S_{H_4} ⊂ I` is enumerated by exhaustive
search inside the bounded coordinate box dictated by the standard
`600`-cell realization; reflection invariance is verified by computing
`r_α(β)` for every pair `(α, β) ∈ S_{H_4}²`.

**Script.** [`scripts/golden_octonions.py`](../scripts/golden_octonions.py)

**Certificates.** [`certificates/icosian_double/icosian_H4_shell_certificate.json`](../certificates/icosian_double/icosian_H4_shell_certificate.json),
[`certificates/icosian_double/icosian_double_root_shell_certificate.json`](../certificates/icosian_double/icosian_double_root_shell_certificate.json)

## (P3) Gram matrix and self-duality

**Content.** Computation of the 4×4 Gram block `G_I` of the polar pairing
`B(x, y) = 2 Re(x ȳ)` on the icosian basis; verification that the entries
agree bit-by-bit with those displayed in equation `eq:gram-icosian` of the
article; computation of `det G_I = φ² = φ + 1`, of the full
`det G_{G_0} = φ⁴ = 2 + 3 φ`, of its `K/Q`-field norm `1`, and of
`G_{G_0}^{-1} ∈ M_8(Z[φ])`.

**Scripts.** [`scripts/g3_dual_discriminant.py`](../scripts/g3_dual_discriminant.py),
[`scripts/g3_e8e8_unimodularity.py`](../scripts/g3_e8e8_unimodularity.py)

**Certificates.** [`certificates/g3/dual_discriminant_certificate.json`](../certificates/g3/dual_discriminant_certificate.json),
[`certificates/g3/e8e8_unimodularity_certificate.json`](../certificates/g3/e8e8_unimodularity_certificate.json)

## (P4) Denominator-two `F_4^8` projective lines

**Content.** Enumeration of all `(4⁸ - 1)/(4 - 1) = 21845` projective lines
of `(Z[φ]/2 Z[φ])^8 ≅ F_4^8`; classification by the three a priori failure
filters

- **(F1)** not mixed (170 lines, = `2 · 85` from the two `(F_4)⁴`-halves);
- **(F2)** mixed but conjugation-stable failure (16320 lines);
- **(F3)** mixed conjugation-stable but integrally unpaired (5355 lines).

The three counts sum to `21845`; survivor count zero.

**Scripts.** [`scripts/g3_search.py`](../scripts/g3_search.py),
[`scripts/g3_candidate_ideals_from_discriminant.py`](../scripts/g3_candidate_ideals_from_discriminant.py)

**Certificates.** [`certificates/g3/denominator2_single_line_search.json`](../certificates/g3/denominator2_single_line_search.json),
[`certificates/g3/denominator2_all_subspaces_no_go.json`](../certificates/g3/denominator2_all_subspaces_no_go.json),
[`certificates/g3/mixed_half_root_search.json`](../certificates/g3/mixed_half_root_search.json)

## (P5) Ramified `√5`-denominator scan

**Content.** Enumeration of all `(5⁸ - 1)/(5 - 1) = 97656` projective lines
of `(Z[φ]/√5)^8 ≅ F_5^8`. Structural argument: the reduced polar form
`B mod √5` is non-degenerate (its discriminant `2 + 3 φ` has nonzero
`K/Q`-field norm modulo `5`), so the integral-pairing filter selects only
the zero line. Survivor count zero.

**Script.** [`scripts/g3_denominator_sqrt5_scan.py`](../scripts/g3_denominator_sqrt5_scan.py)

**Certificates.** [`certificates/g3/denominator_sqrt5_single_line_search.json`](../certificates/g3/denominator_sqrt5_single_line_search.json),
[`certificates/g3/denominator_sqrt5_all_subspaces_no_go.json`](../certificates/g3/denominator_sqrt5_all_subspaces_no_go.json)

## (P6) Trace-polar discriminant tower G3-B'

**Content.** Computation of `G_{0,Z}^# / G_0 ≅ (Z/5 Z)^8`; classification of
the induced quadratic form as split type `O^+(8, 5)` (hyperbolic rank 4);
enumeration of its `19656` isotropic projective lines; verification that
each generates a closure of full dimension under conjugation and
multiplication by the chosen `Z[φ]`-basis of `G_0`. Survivor count zero.

**Scripts.** [`scripts/g3_trace_integral_search.py`](../scripts/g3_trace_integral_search.py),
[`scripts/g3_b_prime_phase2.py`](../scripts/g3_b_prime_phase2.py)

**Certificates.** [`certificates/g3_b_prime/discriminant_quadratic_form.json`](../certificates/g3_b_prime/discriminant_quadratic_form.json),
[`certificates/g3_b_prime/h4_ansatz_search.json`](../certificates/g3_b_prime/h4_ansatz_search.json),
[`certificates/g3_b_prime/h4_ansatz_summary.json`](../certificates/g3_b_prime/h4_ansatz_summary.json),
[`certificates/g3_b_prime/sqrt5_ansatz_search.json`](../certificates/g3_b_prime/sqrt5_ansatz_search.json),
[`certificates/g3_b_prime/sqrt5_ansatz_summary.json`](../certificates/g3_b_prime/sqrt5_ansatz_summary.json),
[`certificates/g3_b_prime/octonion_stable_subspaces.json`](../certificates/g3_b_prime/octonion_stable_subspaces.json)

## Reproducibility

Run from the repository root:

```bash
python verify.py --clean
```

The driver regenerates the certificate tree from scratch and prints a
SHA-256 digest of the final tree. The same digest is stored, per file, in
[`certificates/SHA256SUMS`](../certificates/SHA256SUMS).
