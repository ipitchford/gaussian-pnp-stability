> **Pre-revision record (20 October 2026).** This audit concerns the manuscript before the 27 October 2026 revision. Section numbers, constants (for example the conditioning bound 21, now 19.22) and unqualified "sharp dimension" wording are superseded by the revised gaussian_pnp_stability.md and REVISION_RESPONSE.md wherever they differ.

# Verification record

**Date:** 20 October 2026.  
**Material reviewed:** `gaussian_pnp_stability.md`, `verify_exact.py`,
`verify_circulant.py`, and the resulting rational certificate.  
**Outcome:** no material mathematical gap was found in this review. This is a
separate audit of the assembled manuscript and its computations, not formal
proof-assistant verification or external peer review. Rerunning the supplied
scripts establishes reproducibility; it is not, by itself, an independent proof.

## Exact computations reproduced

Both scripts were executed successfully using Python's standard library.

- `verify_exact.py` reconstructed the four-dimensional matrices from their
  defining formulas and expanded the reduced determinant using rational
  arithmetic. It verified the complete cubic factor and the coefficient
  inequality establishing negativity for **every** $0<t\leq1/2$. It also
  checked the full and reduced determinants at five rational steps, the
  reversed-product determinant identity, all eight coefficients of the
  two-parameter determinant, and the stated mixing-gap and first-loss
  derivatives.
- `verify_circulant.py` verified the positive circulant matrix, its Fourier
  projection identity, the formula for its square, the kernel's detailed
  balance relation, the stationary block weights, and the negative compressed
  coefficient $c=-6287/8364672$. At $t=1/100$,
  $\epsilon=10^{-9}$, it obtained exact negative rational values of
  both $\det(I-P)$ and $\det(I-R)$. The complete fractions are recorded in
  `circulant_certificate.json`.

The twelve-dimensional finite witness is an additional example. It does **not**
lie on the asymptotic path $\epsilon=t^2$, whose instability is established
analytically in the manuscript.

## Analytic arguments checked

The low-dimensional result was checked directly: the three-state positivity
identity makes the rank-one correction of $W$ strictly positive and doubly
stochastic; the correction leaves the quotient action unchanged. Strict
Euclidean contraction of this quotient gives the full claimed interval for
$P$. The spectral-multiplier argument for $R$ is valid also when the Hessian
has a nontrivial nullspace.

For the four-dimensional construction, the audit checked the eigenvalue
decomposition of $A$ and $B$, fixed row-degree ratios and stationary vector,
positive definiteness and conditioning bound for $K_t$, the invariant
antisymmetric direction, and the sign argument converting a negative
determinant at one into a real eigenvalue greater than one. The identity

\[
\det(I-R)=\frac{\det(I-P)}{\det(I+tB)}
\]

was checked from the operator definitions. Its use is legitimate because the
denominator is positive; it does not assert equality of the numerator matrices.

For the circulant construction, the audit checked the block degrees and
stationary weights, the negative two-group compression, the direct-sum
decomposition at zero coupling, and the Schur-complement expansion

\[
\det(I-P)=t^2c\det(I-Q)+O(t^3).
\]

The ten eigenvalues of $Q$ lie in $(0,1)$, so the leading coefficient is
strictly negative. Smooth row normalisation justifies substituting
$\epsilon=t^2$. This proves the sufficiently-small-step statement without
relying on numerical eigenvalues.

The gap-dependent sufficient bound and its first-loss asymptotics were also
checked. A minor notation clarification was passed to the manuscript
maintainer: when the quotient is complexified, its weighted norm minimises over
$c\in\mathbb C$, rather than $c\in\mathbb R$. This does not change the
argument or the bound.

## Gaussian geometry and numerical checks

Both Gaussian representations were verified analytically by computing the
pairwise squared distances. In the four-point construction, the stated bounds
make every square root real. In the twelve-point construction, the orthogonal
private coordinates, shared high-similarity coordinates, and component
separation give exactly the three required off-diagonal affinity levels for
every $0<\epsilon<1$.

Floating-point reconstruction was used only as a supplementary sanity check.
The maximum entry errors reported on rerunning the scripts were below
$4.2\times10^{-17}$ for the sampled four-point kernels and
$7.0\times10^{-18}$ for the twelve-point finite witness. These numerical
checks do not establish the continuum or asymptotic claims; the analytic
arguments and exact rational signs do.

## Scope retained

The sharp dimension claim requires the Hessian to preserve constants up to a
positive scalar. The twelve-dimensional forward matrix is a positive symmetric
periodic convolution; its weights are not assumed Gaussian or decreasing with
distance. The constructions provide Gaussian feature vectors, without claiming
that those vectors arise as overlapping patches of a natural image. The kernel
is selected for the prescribed step and is then fixed during iteration. The
uniform conditioning statements concern the small-coupling regime.

This audit does not establish publication priority, practical prevalence of
divergence, or an editorial quality assessment. Those questions remain distinct from the
mathematical and computational validity checked here.
