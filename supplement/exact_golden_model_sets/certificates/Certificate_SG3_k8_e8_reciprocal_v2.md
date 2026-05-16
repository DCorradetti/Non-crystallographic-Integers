# Certificate: K8/E8 reciprocal distinction v2

## Statement

For the trace-normalized icosian lattice K8,

```text
K8^#/K8 = (Z/5Z)^4.
```

For the standard even unimodular E8 lattice,

```text
E8^#/E8 = 0.
```

Equivalently, the K8 reciprocal basis has fifth denominators, while the E8
reciprocal basis is integral.

## Notation

K8 is the underlying Z-lattice of the icosian order in the basis

```text
1, i, h, g, phi, phi*i, phi*h, phi*g
```

with trace-polar form `Tr_{K/Q}(2*<x,y>_K)`.

## Dependencies

- `3. Certificates/k8_vs_e8/icosian_trace_gram.json`;
- `script_exact_golden_model_sets_v2.py`;
- active run `run_exact_golden_model_sets_20260516_071050`.

## Proof

The certified K8 Gram matrix has determinant `625`.  Its Smith normal form
diagonal is

```text
1, 1, 1, 1, 5, 5, 5, 5.
```

The discriminant quotient of an integral lattice with Gram matrix G is the
cokernel of the embedding determined by G, hence the Smith factors give
`(Z/5Z)^4`.  The computed inverse matrix has denominator lcm `5` and `56`
nonintegral entries.

For the standard E8 Gram matrix, the determinant is `1`, the Smith diagonal
is eight copies of `1`, and the inverse matrix is integral.  Therefore the
dual quotient is trivial.

## Verification of hypotheses

The active JSON comparison file reports:

- K8 determinant: 625;
- K8 Smith diagonal: `[1,1,1,1,5,5,5,5]`;
- K8 dual denominator lcm: 5;
- E8 determinant: 1;
- E8 Smith diagonal: `[1,1,1,1,1,1,1,1]`;
- E8 dual denominator lcm: 1.

## Computational support, if any

Canonical files:

- `reciprocal_module_h4.json`;
- `k8_e8_projection_comparison.json`.

## Failure modes / limitations

The theorem distinguishes integral normalizations and reciprocal modules.  It
does not by itself classify all possible icosian model sets or their atomic
surfaces.

## Public-manuscript status

The proof is included in Article 17 v2.  The certificate supplies the exact
run and JSON trace.

## Changelog

Version 2 adds the Smith normal form diagonal explicitly to the public JSON
comparison and to the manuscript proof.
