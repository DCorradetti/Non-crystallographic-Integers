---
title: Paper ↔ certificate correspondence
---

# Paper ↔ certificate correspondence

This page maps each numerical claim in `NC_Integers_v3.tex` to the JSON
certificate(s) that verify it.

| Paper element                                              | Section / label                        | Certificate                                                                 |
|------------------------------------------------------------|----------------------------------------|-----------------------------------------------------------------------------|
| `Z[φ]` arithmetic baseline (P1)                            | §2, §5                                 | [`field_certificate.json`](../certificates/icosian_double/field_certificate.json) |
| Icosian ring `I` is a `Z[φ]`-order                         | §6.2, Prop. `prop:icosian-h4`          | [`icosian_order_certificate.json`](../certificates/icosian_double/icosian_order_certificate.json) |
| 120-element `H_4` shell, reflection-invariant              | §6.2, Prop. `prop:icosian-h4`          | [`icosian_H4_shell_certificate.json`](../certificates/icosian_double/icosian_H4_shell_certificate.json) |
| `j`, `k` as `Z[φ]`-linear combinations of `e_1,...,e_4`    | Remark `rem:icosian-basis-generates-jk` | [`icosian_order_certificate.json`](../certificates/icosian_double/icosian_order_certificate.json) (basis transformation block) |
| Icosian double `G_0 = I ⊕ Iℓ` is a `Z[φ]`-order            | §6.3, Def. `def:icosian-double`        | [`icosian_double_order_certificate.json`](../certificates/icosian_double/icosian_double_order_certificate.json) |
| 240-element `H_4 ⊕ H_4` shell                               | §6.3, Prop. `prop:weak-golden-octonion` | [`icosian_double_root_shell_certificate.json`](../certificates/icosian_double/icosian_double_root_shell_certificate.json) |
| Nonzero associator `[i,j,ℓ] = 2kℓ`                          | §6.3, eq. `eq:nonzero-associator`      | [`octonion_multiplication_certificate.json`](../certificates/icosian_double/octonion_multiplication_certificate.json) |
| Gram matrix `G_I`, `det G_I = φ²`                          | §6.4, eq. `eq:gram-icosian`            | [`dual_discriminant_certificate.json`](../certificates/g3/dual_discriminant_certificate.json) (`gram_block` field) |
| `det G_{G_0} = φ⁴ = 2+3φ`, field norm = 1                  | §6.4, eq. `eq:gram-g0`                 | [`dual_discriminant_certificate.json`](../certificates/g3/dual_discriminant_certificate.json) (`determinant` and `norm_K_Q` fields) |
| `G_0^# = G_0` (self-duality)                               | §6.4, Prop. `prop:dual-discriminant`   | [`dual_discriminant_certificate.json`](../certificates/g3/dual_discriminant_certificate.json) (`G0_dual_equals_G0`) |
| Z-trace lattice is 5-modular, det = `5⁸`                   | §6.4 final remark; §7.1 remark         | [`e8e8_unimodularity_certificate.json`](../certificates/g3/e8e8_unimodularity_certificate.json) |
| 21845 `F_4`-lines; structural breakdown (F1 + F2 + F3)     | §6.4, Prop. `prop:bounded-g3-no-go`    | [`denominator2_single_line_search.json`](../certificates/g3/denominator2_single_line_search.json), [`denominator2_all_subspaces_no_go.json`](../certificates/g3/denominator2_all_subspaces_no_go.json) |
| 97656 `F_5`-lines; non-degeneracy of `B mod √5`            | §6.4, Prop. `prop:ramified-denominator` | [`denominator_sqrt5_single_line_search.json`](../certificates/g3/denominator_sqrt5_single_line_search.json), [`denominator_sqrt5_all_subspaces_no_go.json`](../certificates/g3/denominator_sqrt5_all_subspaces_no_go.json) |
| 14400 raw pairs, 3600 cosets, all of norm `1/2`            | §6.4, Prop. `prop:mixed-half-root`     | [`mixed_half_root_search.json`](../certificates/g3/mixed_half_root_search.json) |
| `Tr_{K/Q} N(φv) = 3/2` failure for `H_4` mixed half-roots  | §6.4, Prop. `prop:h4-mixed-trace`      | [`h4_ansatz_search.json`](../certificates/g3_b_prime/h4_ansatz_search.json), [`h4_ansatz_summary.json`](../certificates/g3_b_prime/h4_ansatz_summary.json) |
| `(Z/5Z)^8` trace-polar quotient, `O^+(8,5)` type           | §6.4, Prop. `prop:g3-b-tower-no-go`    | [`discriminant_quadratic_form.json`](../certificates/g3_b_prime/discriminant_quadratic_form.json) |
| 19656 isotropic projective lines, full-dim stable closure  | §6.4, Prop. `prop:g3-b-tower-no-go`    | [`octonion_stable_subspaces.json`](../certificates/g3_b_prime/octonion_stable_subspaces.json), [`sqrt5_ansatz_search.json`](../certificates/g3_b_prime/sqrt5_ansatz_search.json) |

## Notes on certificate field conventions

- Counts (`lines_total`, `lines_mixed`, `survivors`, etc.) are recorded as
  Python `int`, exact.
- Algebraic numbers in `Z[φ]` appear as `[a, b]` pairs standing for
  `a + b φ`.
- Algebraic numbers in `K = Q(√5)` appear either as the same pair form, or
  as a rational pair `[a/b, c/d]` standing for the `Z[φ]`-coordinates.
- Field-norm and trace values are stored as plain integers/rationals.
- The Gram matrix `G_I` is stored as a 4×4 nested list with `Z[φ]` entries
  in pair form; the human-readable expression at the bottom of the file
  matches equation `eq:gram-icosian` in the article.

## Hash check

To check that your local certificate tree matches the published one, run

```bash
python verify.py
```

and compare the printed SHA-256 digest with the one listed in the article's
*Computational methods* section (filled in at the camera-ready stage).
