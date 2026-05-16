# Certificate: H2 decagonal benchmark v2

## Statement

The H2 control package constructs a finite cyclotomic model-set benchmark
over the decagonal module and exports exact Gram data, the centered
decagonal strip window, model-set labels, reciprocal labels and finite-patch
kinematic intensities.

## Notation

The source lattice is `Z^5 / diagonal Z`.  Physical columns are
`v_j=(cos(2*pi*j/5), sin(2*pi*j/5))`; internal columns are
`w_j=(cos(4*pi*j/5), sin(4*pi*j/5))`.

The window is the centered zonotope

```text
W_10 = { sum_j t_j w_j : -1/2 <= t_j <= 1/2 }.
```

## Dependencies

- `script_exact_golden_model_sets_v2.py`
- `nc_quasicrystal_core.py`, core version `v2_decagonal_h2_window`
- active run `run_exact_golden_model_sets_20260516_071050`

## Proof

The trigonometric identities for fifth roots put the physical and internal
Gram matrices in `Q(sqrt(5))`.  The window is the convex hull of the 32 sums
`sum_j epsilon_j w_j/2`, with `epsilon_j = +/-1`; the hull has ten vertices,
exported counter-clockwise in `cut_project_basis.json`.

The finite patch is obtained by enumerating integer labels with diagonal sum
zero, applying physical and internal projections, retaining the labels whose
internal projection lies in the decagon, and bounding the physical radius.
The finite intensities use the direct finite sum
`|sum exp(-i k.x)|^2/N^2`.

## Verification of hypotheses

The active run reports:

- H2 patch vertices: 31;
- finite diffraction peaks exported: 25;
- window type: centered regular decagonal zonotope;
- all core checks passed: true.

## Computational support, if any

Canonical files:

- `cut_project_basis.json`;
- `penrose_patch_vertices.json`;
- `reciprocal_module_h2.json`;
- `diffraction_intensities_h2.json`.

## Failure modes / limitations

This is still a finite-window benchmark.  It is not an infinite-volume
diffraction theorem and does not classify Penrose tilings.

## Public-manuscript status

Included in Article 17 v2 as the H2 control computation.

## Changelog

Replaces the v1 circular proxy window with the centered decagonal strip
window.
