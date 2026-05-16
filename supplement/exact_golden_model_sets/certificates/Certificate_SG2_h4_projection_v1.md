# Certificate SG2 - H4 Projection Data v1

## Statement

The H4 package exports the 120 icosian shell coordinates, the projection
convention used for icosahedral star labels, and the twelve standard H3
star vectors over `Z[phi]`.

## Notation

Icosians are represented as quaternions `(a,b,c,d)` over `Q(sqrt(5))`.
The physical projection used in this first implementation is the
imaginary-part map `(a,b,c,d) -> (b,c,d)`, with the scalar coordinate
recorded as a complementary internal coordinate.

## Verification

The run writes:

- `h4_shell_coordinates.json`;
- `projection_matrices_h4.json`;
- `icosahedral_star_vectors.json`.

The exported shell has `120` roots, all with norm one. The star-vector
file exports the twelve vectors
`(0, +/-1, +/-phi)`, `(+/-1, +/-phi, 0)`, and
`(+/-phi, 0, +/-1)`.

## Limitations

This is projection and reciprocal-label infrastructure. It is not yet a
full atomic-surface classification for all icosahedral model sets derived
from the icosian lattice.

## Public-Manuscript Status

Included in Article 17 as the exact H4 coordinate/projection layer.
