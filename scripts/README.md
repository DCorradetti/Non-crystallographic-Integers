# Golden octonions exact checks

This directory contains the first exact-arithmetic implementation required by
`6. Reviews/1/Golden_Octonions_Existence_Roadmap.md`.

Current scope:

- exact arithmetic in `Z[phi]` and `Q(phi)`;
- Cayley-Dickson octonions over `Q(phi)`;
- the icosian ring with basis
  `1, i, (1+i+j+k)/2, (-1+(phi-1)*i-phi*j)/2`;
- the octonionic doubled order `I + I*l`;
- the finite shell `H4 + H4` of 240 norm-one elements;
- exact certificates for closure, conjugation, root-shell checks and a
  nonzero associator.

Run from this directory:

```powershell
python .\golden_octonions.py --out "..\..\3. Certificates\golden_octonions\icosian_double"
```

The current result is a positive, circumscribed construction: the icosian double
is a full-rank `Z[phi]` octonion order with a distinguished non-crystallographic
root shell of type `H4 + H4`. It is not a classification of all possible golden
octonion orders.

## G3 bounded searches

Run:

```powershell
python .\g3_search.py --out "..\..\3. Certificates\golden_octonions\g3"
```

This produces:

- `g1_baseline_invariants.json`
- `denominator2_single_line_search.json`
- `denominator2_all_subspaces_no_go.json`
- `mixed_half_root_search.json`
- `summary.json`

The current G3 slice is negative: no nonzero denominator-2 `F4`-subspace can
pass the strict integral-pairing condition, and the direct mixed half-root
ansatz with `a,b in H4` has norm `1/2`, not an integral `Z[phi]` norm.

## G3 ideal-denominator overorders

Run:

```powershell
python .\g3_dual_discriminant.py --out "..\..\3. Certificates\golden_octonions\g3"
python .\g3_candidate_ideals_from_discriminant.py --cert-dir "..\..\3. Certificates\golden_octonions\g3"
```

This produces:

- `dual_discriminant_certificate.json`
- `candidate_denominator_ideals.json`
- `summary_g3_ideal_search.json`

The norm-polar determinant is `2+3*phi`, a unit of `Z[phi]`. Hence `G0#=G0`
for the norm-polar dual, the candidate denominator ideal list is empty, and
strict overorders of the baseline weak double are excluded.

## G3 strict/relaxed roadmap checks

Run:

```powershell
python .\g3_e8e8_unimodularity.py --out "..\..\3. Certificates\golden_octonions\g3"
python .\g3_denominator_sqrt5_scan.py --out "..\..\3. Certificates\golden_octonions\g3"
```

This produces:

- `e8e8_unimodularity_certificate.json`
- `denominator_sqrt5_single_line_search.json`
- `denominator_sqrt5_all_subspaces_no_go.json`

The trace-polar `Z`-lattice check triggers the normalization-risk branch:
the determinant is `625` on `I` and `5^8` on `G0`, not `1`. The strict
no-go theorem is therefore proved from the `Z[phi]` dual certificate instead
of an `E8 + E8` unimodularity argument. The ramified denominator `(sqrt5)`
has zero projective lines passing the necessary polar-pairing filter.

## G3-B' relaxed H4 ansatz

Run:

```powershell
python .\g3_trace_integral_search.py --out "..\..\3. Certificates\golden_octonions\g3_b_prime"
```

This produces:

- `h4_ansatz_search.json`
- `h4_ansatz_summary.json`

The mixed half-root ansatz has no `Z[phi]`-module survivor under the relaxed
quadratic trace convention. All raw half-roots have `Tr(N(v))=1`, but adjoining
the `Z[phi]`-line also adjoins `phi*v`, with first failure
`Tr(N(phi*v))=3/2`.

## G3-B' relaxed discriminant-tower phase 2

Run:

```powershell
python .\g3_b_prime_phase2.py --out "..\..\3. Certificates\golden_octonions\g3_b_prime"
```

This produces:

- `discriminant_quadratic_form.json`
- `sqrt5_ansatz_search.json`
- `octonion_stable_subspaces.json`
- `sqrt5_ansatz_summary.json`

The trace-polar discriminant quotient is `G0_Z# / G0 ~= (Z/5)^8`. Its
quadratic form is split plus type `O^+(8,5)` and has `19656` isotropic
projective lines. Every isotropic line generates the whole 8-dimensional
quotient under conjugation and left/right multiplication by the `G0` basis,
and the whole quotient is not totally isotropic. Thus the `sqrt5`
trace-discriminant tower has no nonzero octonion-stable isotropic subspace.
