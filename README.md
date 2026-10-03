# Kuramoto finite-step graphon persistence

This repository develops the graphon side of the six-block Kuramoto program.
It is separate from the direct finite blow-up proof in
[`kuramoto-block-limit`](https://github.com/ssiddhantsood/kuramoto-block-limit)
because the two projects have different hard steps:

- this repository studies infinite-dimensional step-graphon operators and
  persistence under graph approximation;
- the direct repository studies exact equitable finite blow-ups and their
  transverse spectra.

## First result

For a finite-step symmetric graphon and a phase function constant on its
cells, the Kuramoto linearization splits exactly into:

1. a finite-dimensional mass-orthonormal quotient operator; and
2. cell-zero-mean modes on which the operator is multiplication by `-q_i`.

Consequently, linear stability modulo rotation is equivalent to positivity of
the nonrotation quotient Hessian eigenvalues and positivity of every `q_i`.
The statement and proof are in
[`docs/STEP_GRAPHON_SPECTRUM.md`](docs/STEP_GRAPHON_SPECTRUM.md).

For the current near-boundary six-block candidate, run:

```bash
python3 scripts/check_candidate.py
python3 -m unittest discover -s tests -v
```

The calculation finds two quotient modes near `1e-12`, while the three
distinct within-cell quantities remain comfortably positive.

## Persistence target

The next theorem adapts Bramburger--Holzer--Williams (arXiv:2402.09276) from
continuous graphons/equilibria to a fixed finite partition with jumps. See
[`docs/PERSISTENCE_EXTENSION.md`](docs/PERSISTENCE_EXTENSION.md).

