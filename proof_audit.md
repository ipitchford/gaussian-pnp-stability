> **Pre-revision record (20 October 2026).** This audit concerns the manuscript before the 27 October 2026 revision. Section numbers, constants (for example the conditioning bound 21, now 19.22) and unqualified "sharp dimension" wording are superseded by the revised gaussian_pnp_stability.md and REVISION_RESPONSE.md wherever they differ.

---
title: "Proof audit and human-readable outline"
subtitle: "Track I: informal mathematics with exact reproducibility checks"
date: "20 October 2026"
lang: en-GB
---

## Audit conclusion

The main manuscript's three theorems and the first-loss asymptotics pass the analytical and computational checks recorded here. No unresolved substantive mathematical gap was identified. The proof is unusually explicit about its quantifiers, its spectral conclusion, and the limits of its imaging interpretation.

This is an internal collaborative audit, including separate rederivations of key identities. It is not formal proof-assistant verification or independent external peer review. Rerunning the supplied script establishes reproducibility; it is not by itself an independent proof. Novelty is assessed in the separate novelty report.

The scope here is *gaussian_pnp_stability.md*, its Appendix A, and the two main verifiers. The later three-state reentrance supplement is a separate output and is not included in the claim ledger or rubric below. The separate numerical stability diagram is explanatory; no theorem depends on it.

## 1. Statement

The two matrices are
\[
P_t=W(I-tB),\qquad
R_t=I-W+(I+tB)^{-1}(2W-I),
\]
where \(W=D^{-1}K\), \(D=\operatorname{diag}(Ke)\).

The central classification says that a connected symmetric unit-diagonal kernel with entries in \([0,1]\), together with a symmetric positive semidefinite Hessian satisfying \(Be=b e\), \(b>0\), gives stability in dimensions at most three throughout \(0<t<2/\rho(B)\). In dimension four this can fail for genuine Gaussian kernels, even when the strictly positive averaging forward matrix, stationary vector and conditioning bounds are fixed.

A complementary theorem guarantees stability when the step is small relative to the mixing gap. The counterexample family shows that this dependence has the correct asymptotic order. A separate twelve-state construction gives the failure with a positive circulant forward matrix.

The phrase “for every prescribed step” means that an instance is selected for that step. A single kernel remains fixed during its iteration. The result does not assert instability at every small step for one fixed positive kernel.

## 2. Idea

The geometric difficulty is that the data update and the denoiser contract in different geometries. The Hessian is positive in the Euclidean inner product. The denoiser averages according to its stationary measure. Near a disconnected affinity graph, vectors that are constant inside each component form a slow subspace.

Compressing the Hessian onto component levels uses one map to assign constants and another to take weighted averages. Those maps are not transposes. The effective Hessian can therefore have a negative eigenvalue even though the original Hessian is positive definite. Mixing damps a component contrast at a rate proportional to the graph gap; the negative compressed coefficient makes the data update amplify it at a rate proportional to the step.

The low-dimensional theorem identifies an obstruction to this mechanism. With at most three vertices, a simple correction of the row-normalised kernel is strictly positive and doubly stochastic, while retaining the same action modulo constants. That common Euclidean quotient contraction makes both algorithms stable for the full usual step interval.

## 3. Proof outline

| Step and location | Mathematical move | Classification |
|---|---|---|
| Separate constants, §3, equation (4) | Use an orthonormal basis beginning with \(e/\sqrt n\). The denoiser is upper triangular and the Hessian is block diagonal; control the scalar constant mode separately. | Standard structure, clearly stated |
| Prove three-state positivity, §3.1, equations (5)–(8) | The correction \(S=W+e(e^{\mathsf T}-e^{\mathsf T}W)/3\) has positive entries by an explicit rational identity. Strict Jensen contraction gives \(\|Q\|_2<1\) on the quotient. | The main point of the dimension theorem |
| Transfer quotient contraction to both algorithms, §3.2, equation (9) | A product norm handles ISTA. An eigenvalue-resolvent argument excludes every \(|\lambda|\geq1\) for ADMM. | Short reusable argument |
| Construct a valid four-state family, §4, equations (10)–(16) | A rank-one spectral decomposition fixes \(A,B\). Equal degree vectors of \(K_0,H\) fix \(\pi\). Elementary perturbation bounds control positive definiteness and conditioning. | Explicit construction plus routine checks |
| Certify instability on the whole interval, §5, equations (17)–(22) | Exchange symmetry removes one positive determinant factor. Exact expansion leaves a cubic with a uniformly negative upper bound. A determinant identity transfers the sign to ADMM. | The main point of the counterexample certificate |
| Verify genuine Gaussian geometry, §6, equations (23)–(25) | Construct four feature vectors and compute all six squared distances; analytic inequalities make the square roots real. | Explicit geometric bridge |
| Relate stability to the gap, §7, equations (26)–(36) | A stationary-weight quotient norm gives a sufficient step bound. A two-parameter determinant and a simple eigenvalue branch give the first crossing and its ratio to the gap. | Standard perturbation tools, new explicit family constants |
| Impose convolution, §8, equations (37)–(42) | A rational Fourier projection gives a positive circulant \(A_c\). Weighted two-block compression has a negative eigenvalue. A Schur complement along \(\epsilon=t^2\) gives the sign, and an orthogonal feature construction gives exact Gaussian affinities. | Additional structured construction |

The large integer coefficients are not the proof's conceptual bottleneck. They certify a small matrix identity. The essential choices are the fixed stationary weights, the negative oblique compression, and the positivity identity that rules out the mechanism in dimension three.

## 4. Components and dependencies

### What is supplied explicitly here

The four-state parameter choice, its exact determinant certificate, the feature coordinates, the special three-state positivity argument, the matching-order first-loss family, and the twelve-state circulant construction are all given in full. Their status as original literature contributions remains bounded by the separate novelty search.

### What is standard

The proof uses the spectral theorem, Perron–Frobenius theory for positive matrices, Jensen's inequality, elementary operator-norm comparison, Sylvester's determinant identity, the intermediate value theorem, a Schur complement, and the implicit function theorem for a simple eigenvalue. It does not require a new general theorem about oblique projections or Markov aggregation.

### Dependency order

The low-dimensional proof is independent of the arithmetic construction. The four-state theorem depends on §§4–6. The gap-order obstruction depends on the same family and Appendix A's two-parameter quotient. The circulant theorem has its own construction and limiting proof, sharing only the general determinant bridge for ADMM. This makes failures or revisions locally identifiable.

## 5. Prior work relevant to interpreting the proof

The broad conjecture is stated in Banerjee and Chaudhury (2025), whose fixed-instance local theorem is consistent with this work. Patel's posted example already disproves the unrestricted interval and has arbitrarily small thresholds. Gavaskar et al. (2021, Appendix C) already gives actual kernel ADMM divergence with a signed measurement and an unequal-diagonal kernel.

Sinha et al. (2025) supplies the unit-diagonal kernel setting and distinguishes standard from scaled algorithms; its Gaussian norm obstruction does not itself prove spectral instability. Gavaskar and Chaudhury (2020) and Liu et al. (2021) contain positive convergence results with additional conditions. Liu et al.'s complete theorem could not be retrieved during the bounded search.

The proof does not depend on treating any of these sources as absent. Its theorem statements specify the restrictions that distinguish its contribution. Full citations, inspected sections, exact translations, and access limitations are in *novelty_report.md* and the eleven-part *academic_review.md*.

## 6. Verification record

### Environment and execution

The two supplied scripts were inspected before execution. They use Python's standard library for exact rational arithmetic. The final copies were rerun on 20 October 2026 with Python **3.12.14**. Both exited successfully with empty standard error. Complete commands, working directory, elapsed times and output are retained in *proof_audit_verification.txt*.

| Script | Exit | Elapsed time | SHA-256 of the executed file |
|---|---:|---:|---|
| verify_exact.py | 0 | 0.031423 s | f96d218a0cd06b90c8f5b8ec44eaef7a91dc191137fe142e4e42e08366df5f7f |
| verify_circulant.py | 0 | 0.054699 s | 2cd488f34f589d8f48b69073f8cc7c2b8c3beb0b88bf5ffbbce0b70b47ce3f34 |

The reported timings are observations from this execution, not performance claims.

### What the scripts establish

The four-state verifier reconstructs the matrices, expands the reduced determinant exactly, and verifies the coefficient inequality making its cubic negative throughout \(0<t\leq1/2\). It also checks selected rational instances, the determinant bridge, all eight nonzero coefficients of the bivariate determinant, and the gap and first-loss derivatives.

The circulant verifier checks the Fourier projection identities, positive forward matrix, detailed balance, component weights, negative compressed coefficient, and exact negative determinant signs for the additional finite witness \(t=1/100,\epsilon=10^{-9}\). It writes the complete rational values to *circulant_certificate.json*.

Floating-point reconstruction of Gaussian coordinates is supplementary. The maximum sampled entry error was below \(4.2\times10^{-17}\) for the four-state family and \(7.0\times10^{-18}\) for the circulant witness. These samples do not establish the continuum or asymptotic conclusions.

### Separate analytical and algebraic checks

A second calculation reconstructed all eight bivariate coefficients from the original matrices using rational arithmetic, separately from the final verifier. It confirmed
\[
g'_0=\frac{31}{26},\quad
\partial_t\lambda(0,0)=\frac{133}{812500},\quad
T'_0=\frac{968750}{133},\quad
\lim_{\epsilon\to0}\frac{T_\epsilon}{g_\epsilon}=\frac{812500}{133}.
\]

The reviewer independently expanded all six four-point Gaussian distances and the three distance types in the hierarchical construction. The review also checked the three-state positivity identity, the ADMM resolvent bound, the norm-comparison constants, and the circulant Schur-complement argument. These are the analytic parts that a successful script run alone would not establish.

## 7. Claim ledger

“Verified” in this ledger means checked within this session, not certified in a proof assistant.

| Claim | Location | Evidence | Status and limitation |
|---|---|---|---|
| Constant and quotient modes are separated correctly | §3, (4) | Direct block calculation using \(We=e\), \(Be=b e\), and symmetry of \(B\) | Verified; no false assumption that \(e^\perp\) is invariant under \(W\) |
| Three-state correction is strictly positive and doubly stochastic | §3.1, (5)–(7) | Algebraic expansion; direct-edge or two-edge-path argument; \(ab_0\leq1\) for the diagonal bound | Verified under connectedness and entries at most one |
| Strict Euclidean quotient contraction | §3.1, (8) | Strict Jensen inequality on nonconstant vectors and compactness | Verified; not a general statement about arbitrary stochastic matrices |
| Full step interval for both operators in \(n\leq3\) | §3.2 | Product norm for \(P\); spectral multiplier for \(R\); scalar modes checked separately | Verified, including singular \(B\) with \(b>0\) |
| Fixed positive stochastic \(A,B\) and their spectra | §4.1, (10)–(12) | Orthogonal decomposition into \(e,v,\{e,v\}^{\perp}\) | Verified exactly |
| Fixed \(\pi\), positive spectrum of \(W_t\), \(\kappa_2(K_t)<21\) | §4.2, (13)–(16) | Equal degree vectors, similarity to a symmetric positive definite matrix, elementary spectral bounds | Verified on the declared interval |
| Spectral instability for every prescribed \(0<t\leq1/2\) | §5, (17)–(21) | Positive antisymmetric factor and exact uniformly negative cubic factor | Verified continuum conclusion, not inferred from sampled steps |
| Transfer of instability to \(R_t\) | §5, (22) | Direct multiplication, Sylvester identity, positive determinant denominator | Verified; no commutativity assumed |
| Exact Gaussian realisation in four states | §6, (23)–(25) | Analytic radicand bounds and all six distances | Verified; feature diameter grows as coupling vanishes |
| Sufficient gap bounds | §7, (26)–(30) | Self-adjoint stationary-weight quotient and Euclidean-to-weighted norm comparison | Verified; constants conservative |
| First-loss asymptotics and sharp order | §7.1, (31)–(35); Appendix A.2 | Separate exact coefficient check, simple eigenvalue branch, implicit function theorem, remaining modes separated | Verified for sufficiently small positive \(\epsilon\); not a sharp universal coefficient |
| Negative circulant compression and bounded conditioning | §8.1–8.3, (37)–(40) | Exact Fourier and degree calculations; \(c=-6287/8364672\) | Verified; the stationary vector varies in this additional example |
| Circulant small-step instability | §8.4, (41) | Smooth normalisation and Schur complement with positive \(\det(I-Q)\) | Verified asymptotically along \(\epsilon=t^2\); the finite witness is a different parameter pair |
| Hierarchical Gaussian realisation | §8.5, (42) | Fifteen mutually orthonormal directions, private/shared coordinates, component separation | Verified for all \(0<\epsilon<1\); no minimal feature dimension claimed |
| Literature priority | §§1.2, 9; separate reports | Bounded source search and reviewer judgement | Provisional; not a mathematical verification claim |

## 8. Adversarial checks and resolved points

**Norm failure was not substituted for spectral divergence.** The proof establishes a negative characteristic polynomial at one and hence a real eigenvalue exceeding one. It does not infer instability from a large singular value.

**The two geometries were kept separate.** Reversibility gives a weighted self-adjoint operator; it does not make the original row-normalised matrix Euclidean symmetric. The gap theorem explicitly pays the norm-comparison factor \(\sqrt{\kappa_\pi}\).

**The ADMM estimate was checked for complex eigenvectors.** In §3.2 the denominator obeys
\[
|\lambda(1+t\beta)-t\beta|
\geq|\lambda|(1+t\beta)-t\beta\geq1
\]
when \(|\lambda|\geq1\). The numerator satisfies \(|1-t\beta|\leq1\) for \(0<t<2/L\). No firm nonexpansiveness assumption on the quotient is missing. The weighted quotient norm in Theorem 2 now explicitly minimises over complex constants when its vectors are complex.

**The first local crossing is also the first loss for small coupling.** At the limiting point the slow quotient eigenvalue is simple, its step derivative is positive, and all other quotient eigenvalues stay uniformly inside the unit disk. The constant modes are stable for the small positive steps considered. This rules out an earlier nearby crossing by another mode.

**The circulant compression was not inferred by adding arbitrary eigenvalues.** Each relevant compressed matrix preserves the constant vector. On the one-dimensional quotient, the nonconstant coefficients add. The manuscript additionally explains why the compressed all-ones term contributes zero.

**The Schur complement has the stated order.** Its upper block is \(tFB_cU+O(t^2)\), its off-diagonal product contributes only \(O(t^2)\), and its lower block is invertible near zero. The resulting two-dimensional determinant has leading term \(t^2c\), multiplied by \(\det(I-Q)>0\).

**All feature directions are explicitly orthonormal.** This closes a possible grammatical ambiguity in the hierarchical Gaussian construction. It is an exposition clarification; the intended distance calculation was already correct.

**Uniform conditioning is not uniform geometric connectivity.** The smallest cross-component affinity tends to zero and the four-point feature diameter becomes unbounded. The affinity matrices can remain well conditioned while the row-normalised chains mix arbitrarily slowly. The final manuscript states this limitation. No uniform feature-diameter bound is claimed.

**The construction does not prove routine practical divergence.** The feature vectors are not shown to be overlapping image patches, the circulant averaging kernel is not a conventional monotone Gaussian blur, and the spectral excess can be small. These affect application and significance claims, not the exact matrix theorem.

## 9. Fifteen-submetric structure and readability rubric

This is the proof-digest Track I rubric, scored on a **1–5 scale**, allowing half points. For informal mathematics, “tactic transparency” means transparency of the reasoning steps; “broad tactic control” means that routine algebra and computational certificates have a declared scope; and “mathlib style” is interpreted as consistent mathematical notation and exposition. No Lean-specific diagnosis is applied.

The scores concern presentation and proof organisation. They are not correctness probabilities or novelty scores.

| Dimension | Submetric | Score | Evidence |
|---|---|---:|---|
| Structure | Main theorem slimness | 4.0 | Theorem 1 states the central dichotomy and quantifiers, although its several fixed-parameter properties make it long. |
| Structure | Complexity distribution | 4.5 | The positivity idea appears before the construction; large exact matrices and coefficients are isolated in §5 and Appendix A. |
| Structure | Dependency clarity | 4.5 | The reading guide and modular §§3–8 distinguish the dimension theorem, certificate, gap law and circulant extension. |
| Signature | Statement naturalness | 4.0 | Constant preservation is well motivated; the long list of frozen constants is necessary for the obstruction but somewhat specialised. |
| Signature | Binder economy | 4.5 | Assumptions and quantifiers are explicit, and the main normalisation \(\rho(B)=1\) avoids redundant scale parameters. |
| Signature | Generality | 4.5 | The low-dimensional theorem does not unnecessarily require \(K\succeq0\) or entrywise nonnegative \(B\); the \(Be=b e\) condition is retained. |
| Step transparency | Reasoning transparency | 4.5 | §2 explains the slow-mode mechanism; determinant signs and spectral consequences are connected in prose. |
| Step transparency | Explicit lemma use | 4.0 | Standard ingredients are named, and the ADMM multiplier is derived. The perturbation-to-first-loss step remains relatively compressed. |
| Step transparency | Control of routine algebra and computation | 4.5 | The exact polynomial is exposed, the sign proof is short, and numerical checks are labelled supplementary. |
| Reuse | Helper usefulness | 4.5 | The quotient-contraction lemma and determinant bridge serve multiple results rather than one numerical instance. |
| Reuse | Reuse potential | 4.0 | The oblique-compression and gap arguments are reusable; the particular rational construction remains specialised. |
| Reuse | No unused helpers | 5.0 | Each displayed intermediate matrix or identity has a clear role in a theorem or verification certificate. |
| Readability | Proof readability | 4.5 | Context, a reading guide, equation locators and explicit scope statements make the argument navigable. |
| Readability | Mathematical style | 4.5 | Norms, spectra, real/complex quotients and parameter conventions are consistently explained after the minor clarification. |
| Readability | Maintainability | 4.5 | Small standard-library verifiers, exact coefficients, a claim ledger and separate numerical visuals support later revision. |

Dimension means:

| Dimension | Mean |
|---|---:|
| Structure | 4.33 |
| Signature | 4.33 |
| Step transparency | 4.33 |
| Reuse | 4.50 |
| Readability | 4.50 |
| **Overall mean of the five dimensions** | **4.40 / 5** |

The structure score is strongest where the manuscript separates the conceptual slow-mode mechanism from arithmetic. Its statement score reflects deliberate generality in the positive theorem, balanced against a necessarily specialised explicit example. Step transparency is high because the determinant and Gaussian-distance arguments can be checked locally; the first-loss paragraph is the densest remaining step. Reuse is strongest for the common quotient and determinant bridge. Readability benefits from the clear quantifiers and limitations; a specialist audience is still assumed for Schur complements and simple-eigenvalue perturbation.

## 10. Suggested final presentation choices

No new theorem-sized repair is required. Preserve the explanation in §2 before the detailed matrices, retain the explicit distinction between an instance and a family, and keep the short sign bound visible in the main text. Readers should be able to understand why the example works before checking its exact integers.

For a journal version, a single named lemma for “Euclidean quotient contraction implies stability of both operators” would make §3.2 easier to cite. A short separate lemma collecting the negative-compression asymptotic could similarly serve both the intuition and the circulant proof. These are optional organisational improvements; adding many helper lemmas would make this compact paper harder to read.

This audit establishes that the mathematical package is concrete and reproducible enough for editorial assessment; it does not predetermine its outcome.

