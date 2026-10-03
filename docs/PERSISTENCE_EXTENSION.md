# Partition-adapted persistence program

## Relationship to arXiv:2402.09276

Bramburger--Holzer--Williams prove persistence and stability transfer when the
limiting equilibrium is continuous, the graphon satisfies a row-continuity
condition, and the linearization is invertible. Their Section 6.4 sketches a
piecewise-continuous function space for bipartite graphons and states that the
same idea should extend to graphons with finitely many jumps, including block
models.

Our finite-step graphon and class-constant phase violate the global continuity
hypotheses at cell boundaries. Moreover, the conjectured optimal point has two
additional zero modes and therefore violates invertibility even after removing
rotation.

## Proposed theorem

Fix a partition `P={I_1,...,I_k}` and define

```text
X_P = direct-sum_i C(I_i)
```

with the maximum of the cellwise sup norms. Work on its mean-zero subspace to
remove rotation. Assume:

1. `W` is bounded and continuous on every rectangle `I_i x I_j`;
2. graph approximants converge in cut norm and their degree functions converge
   uniformly on every cell;
3. a piecewise-continuous equilibrium `u*` exists;
4. the restricted linearization is invertible; and
5. its spectrum is bounded strictly inside the left half-plane.

Then sufficiently large graph approximants possess nearby equilibria whose
restricted linearizations are also stable.

## Proof ledger

The proof should be organized as an audit of the published argument:

| Component | Required change |
|---|---|
| Banach space | Replace `C[0,1]` with `X_P` cellwise |
| Multiplication operator | Require a positive cellwise lower bound for `q_i` |
| Compact integral term | Prove compactness rectangle by rectangle |
| Inverse | Work on the mean-zero subspace |
| Graphon estimates | Make constants uniform over finitely many rectangles |
| Newton map | Reuse the second-iterate contraction argument |
| Stability transfer | Apply spectral separation on the restricted space |

## Immediate next lemma

Prove that the linearization on `X_P` is a Fredholm operator of index zero when
all cellwise multiplication coefficients are bounded away from zero. For an
exact step graphon this follows transparently from the spectral decomposition
in `STEP_GRAPHON_SPECTRUM.md`; the cellwise-continuous case is a compact
perturbation of an invertible multiplication operator.

## Boundary behavior

The persistence theorem applies to strict interior points. To analyze the
limiting optimizer itself, introduce the minimum-degree ratio as a parameter
and perform a finite-dimensional Lyapunov--Schmidt reduction for the two weak
modes. That is a separate theorem, not a consequence of ordinary persistence.

