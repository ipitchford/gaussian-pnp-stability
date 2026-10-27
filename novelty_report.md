> **Pre-revision record (20 October 2026).** This audit concerns the manuscript before the 27 October 2026 revision. Section numbers, constants (for example the conditioning bound 21, now 19.22) and unqualified "sharp dimension" wording are superseded by the revised gaussian_pnp_stability.md and REVISION_RESPONSE.md wherever they differ.

---
title: "Formula-first novelty audit"
subtitle: "Gaussian-kernel plug-and-play stability, sharp dimension, and gap dependence"
date: "20 October 2026"
lang: en-GB
---

## Decision and scope

**Outcome for the restricted main contribution: BOUNDED UNCERTAINTY.** The bounded search located no verified prior statement subsuming the combination of the sharp dimension threshold, genuine Gaussian affinities, fixed positive averaging forward operator, fixed stationary distribution, bounded affinity conditioning, and vanishing first-loss threshold. This is a defensible candidate contribution; it is not an established priority claim.

**Outcome for the unrestricted disproof: COLLISION.** Patel's posted two-state family already disproves the broad Banerjee–Chaudhury interval and gives arbitrarily small failure thresholds. Gavaskar, Athalye and Chaudhury (2021, Appendix C) also already establish spectral divergence of standard kernel PnP-ADMM with a signed measurement vector. Neither result is claimed as new here.

This report audits the final manuscript *A sharp dimension threshold for Gaussian-kernel plug-and-play stability*, rather than earlier exploratory examples. It separates object identity, mathematical verification, and literature novelty. Exact arithmetic supports the theorem; a failed keyword search cannot establish its novelty.

## 1. Frozen claims and quantifiers

Let \(e\) be the all-ones vector, \(K\) a symmetric affinity matrix, \(D=\operatorname{diag}(Ke)\), \(W=D^{-1}K\), and \(B=A^{\mathsf T}A\). The two unscaled error operators are

\[
P_t=W(I-tB),\qquad
R_t=I-W+(I+tB)^{-1}(2W-I).
\]

The claims under audit are:

1. **Sharp dimension in the constant-preserving class.** If \(n\leq3\), \(K\) is nonnegative, symmetric and irreducible, \(K_{ii}=1\), \(K_{ij}\leq1\), and \(B\) is symmetric positive semidefinite with \(Be=b e\), \(b>0\), then \(P_t\) and \(R_t\) are Schur stable for \(0<t<2/\rho(B)\). Positive semidefiniteness of \(K\) is not needed for this positive theorem. Failure becomes possible in dimension four even with genuine Gaussian affinities and strictly positive \(A\).
2. **A four-state family with fixed parameters.** There exist fixed strictly positive symmetric stochastic positive definite \(A\), fixed \(B=A^2\), and fixed \(\pi=(10,10,6,5)/31\), such that for every prescribed \(0<t\leq1/2\) a genuine Gaussian kernel \(K_t\) has \(\kappa_2(K_t)<21\), while both \(P_t\) and \(R_t\) have a real eigenvalue exceeding one. Here \(\rho(B)=1\), \(\kappa_2(B)=2500\), and the degree ratio is exactly two.
3. **Gap dependence with matching asymptotic order.** For reversible primitive \(W\) with nonnegative spectrum, gap \(g=1-\lambda_2(W)\), and stationary ratio \(\kappa_\pi\), the sufficient step bound \(t<g/\sqrt{\kappa_\pi}\) holds when \(B\) is symmetric positive semidefinite and stochastic. The four-state family has first-loss threshold \(T_\epsilon/g_\epsilon\to812500/133\). This establishes the order in \(g\), not an optimal universal constant.
4. **An additional convolution example.** A fixed strictly positive symmetric circulant stochastic positive definite \(12\times12\) forward matrix, with \(\kappa_2(A)=12\), gives spectral instability of both operators along a genuine Gaussian family \(K_\epsilon\), \(\epsilon=t^2\), for all sufficiently small positive \(t\). Its stationary ratio is bounded by \(50/13\), and \(\kappa_2(K_\epsilon)<413\) for \(\epsilon\leq1/100\).

The quantifier in claim 2 is \(\exists A,\pi\ \forall t\ \exists K_t\). The kernel is fixed during an iteration, but varies between problem instances. The claim is not \(\exists K\ \forall\) sufficiently small \(t\). The positive per-instance local theorem in Banerjee–Chaudhury remains valid.

Excluded claims include first general kernel divergence; the first unrestricted counterexample; dimension four minimality without the constant-preserving Hessian assumption; a sharp numerical universal step constant; a Gaussian or monotone convolution point-spread function; or features realised as overlapping patches of a natural image.

## 2. Normalisation and object-identity table

| Object or convention | Translation into this manuscript | Consequence for comparison |
|---|---|---|
| ISTA step \(\gamma\), least-squares loss \(\frac12\|Ax-y\|^2\) | \(t=\gamma\), \(B=A^{\mathsf T}A\) | The error matrix is exactly \(W(I-tB)\). |
| ADMM penalty \(\rho_{\rm pen}\) in Gavaskar et al. (2021) | \(t=1/\rho_{\rm pen}\) | A smaller \(t\) is a larger penalty in that convention. |
| Lifted ADMM matrix in Gavaskar et al. | If \(\mathcal U=[W;I-W]\), \(\mathcal V=[C^{-1},I-C^{-1}]\), \(C=I+tB\), then \(\mathcal V\mathcal U=I-W+C^{-1}(2W-I)\) | The nonzero eigenvalues of their lifted matrix \(\mathcal U\mathcal V\) coincide with those of the present \(R_t\). This is an actual spectral collision for an unrestricted claim. |
| Symmetric kernel normalisation | \(D^{1/2}WD^{-1/2}=D^{-1/2}KD^{-1/2}\) | Similarity preserves spectrum but changes the Euclidean geometry of the data operator. |
| Stationary measure | \(\pi_i=D_{ii}/\sum_jD_{jj}\) | Multiplying all degrees by one common scalar leaves \(\pi\) unchanged. |
| Changing units of the Hessian | Replace \(B\) by \(B/L\) and \(t\) by \(tL\), \(L=\rho(B)\) | Statements with \(L=1\) must be compared after rescaling. |
| Constants versus the Euclidean complement | \(W\) is upper triangular in a basis \((e/\sqrt n,U)\); its quotient is \(U^{\mathsf T}WU\) | The complement need not be invariant. Confusing an orthogonal compression with a restriction would invalidate a proof. |
| Scaled PnP algorithms | Their gradient/proximal step is changed to use a denoiser-dependent metric | Their convergence theorems concern different operators, not the unscaled \(P_t,R_t\). |
| Norm failure versus spectral instability | \(\|P\|>1\) in a specified norm does not imply \(\rho(P)>1\) | Sinha et al.'s norm counterexample does not establish the spectral conclusion proved here. |

Allowed comparisons include simultaneous vertex permutations, positive common rescaling of \(K\), changes of feature coordinates preserving pairwise distances, and the stated \(B,t\) normalisation. An arbitrary diagonal similarity of the entire iteration is useful for spectra but does not preserve all the claimed Euclidean forward-model restrictions. Replacing \(K\) by a correlation-normalised matrix is also not an identity: it generally changes \(W\).

## 3. Formula and small-instance fingerprints

### 3.1 Four-state construction

The central exact matrix fingerprint is

\[
A=\frac1{1000}
\begin{pmatrix}
329&309&109&253\\
309&329&109&253\\
109&109&554&228\\
253&253&228&266
\end{pmatrix}
=\frac{20I+245ee^{\mathsf T}+vv^{\mathsf T}}{1000},
\quad v=(-8,-8,17,-1)^{\mathsf T}.
\]

Its spectrum is \(\{1,219/500,1/50,1/50\}\). The Hessian spectrum is
\(\{1,47961/250000,1/2500,1/2500\}\).

\[
K_0=\begin{pmatrix}
1&9/10&1/10&0\\
9/10&1&1/10&0\\
1/10&1/10&1&0\\
0&0&0&1
\end{pmatrix},\qquad
H=\begin{pmatrix}
0&3/2&3/10&1/5\\
3/2&0&3/10&1/5\\
3/10&3/10&0&3/5\\
1/5&1/5&3/5&0
\end{pmatrix}.
\]

For \(\epsilon=t/100000\), \(K_t=K_0+\epsilon H\), and the degree vector is
\((1+\epsilon)(2,2,6/5,1)^{\mathsf T}\). The exact choices \(t=1/100,1/4,1/2\) give \(\epsilon=1/10000000,1/400000,1/200000\), respectively. These are three rational instances of the same frozen object. The exact verifiers and the accompanying machine fingerprints permit comparison of spectra and characteristic polynomials without depending on labels or rounded eigenvalues.

The reduced determinant factor is \(t^2(1+t/100000)f(t)\), where the coefficient fingerprint, in increasing powers of \(t\), is

\[
f=\left(
-\frac{1973}{25000000},
\frac{74272809}{500000000000},
-\frac{236896251}{625000000000000000},
-\frac{4748139}{1562500000000000000000}
\right).
\]

Further structural fingerprints are the negative slow compression
\(-133/812500\), the gap derivative \(31/26\), the first-loss derivative
\(968750/133\), and the threshold-to-gap limit \(812500/133\).

### 3.2 Low-dimensional identity

For
\(K=\left(\begin{smallmatrix}1&a&b_0\\a&1&c\\b_0&c&1\end{smallmatrix}\right)\),
\(d_1=1+a+b_0,d_2=1+a+c,d_3=1+b_0+c\), the central identity is

\[
S=W+\frac13 e(e^{\mathsf T}-e^{\mathsf T}W),\qquad
3S_{21}
=\frac{a(1+c)+b_0c}{d_1d_3}+\frac{2a}{d_2}.
\]

The general correction is elementary and is not itself claimed as new. The feature requiring a prior-art check is its strict entrywise positivity for connected three-state unit-diagonal kernels with entries at most one, and the resulting full-interval stability classification.

### 3.3 Circulant construction

The first row of \(A_c\) is

\[
(33,11,1,11,21,11,1,11,21,11,1,11)/144,
\]

with spectrum \(\{1,(1/2)^{[2]},(1/12)^{[9]}\}\). The two-block partition is
\(G=\{0,1,2,4,6\}\), \(H=G^c\), with high-affinity groups
\(\{0,2,4,6\}\) and \(\{5,8,9,10\}\). The Fourier compression coefficient is
\(-12875/406616\); the Hessian compression is \(-6287/8364672\).
The exact finite witness uses \(t=1/100,\epsilon=10^{-9}\), separately from the asymptotic path \(\epsilon=t^2\).

These are finite matrix and polynomial objects. No natural integer counting sequence or generating function has been identified. No OEIS match or exhaustive sequence-database exclusion is claimed.

## 4. Alias map used for the search

| Literature | Aliases and mathematical object |
|---|---|
| Optimisation and imaging | Plug-and-play ISTA, forward-backward splitting, standard/unscaled PnP-ADMM, linear kernel denoiser, NLM, row-normalised Gaussian affinity, convergence counterexample |
| Matrix analysis | Schur stability, products of positive semidefinite contractions in different metrics, diagonally symmetrisable stochastic matrices, rank-one stochastic correction, similarity to doubly stochastic matrices |
| Probability and networks | Reversible Markov kernel, almost reducible chain, two-component aggregation, stationary-weight quotient, slow mixing, consensus disagreement norm |
| Model reduction | Oblique compression, Petrov–Galerkin projection, loss of positivity under nonorthogonal compression, stability preservation |
| Geometry and kernel methods | Gaussian distance matrix, unit-diagonal positive definite kernel, squared Euclidean distances, hierarchical Gaussian feature embedding |
| Convolution | Circulant positive averaging matrix, periodic blur, translation-invariant forward map, Fourier projection |

The search treats these as possible representations of the same mechanism, not as separate claims of novelty. Standard oblique projection instability, norm comparison, Perron–Frobenius theory, and the implicit function theorem receive no novelty credit.

## 5. Retained query log

All searches below were executed on **20 October 2026**. The log records exact retained strings used in this session and their consequences. It is not presented as an exhaustive index of all exploratory browsing. Engine 1 and Engine 2 denote the two available web-search backends. Relevant primary sources were opened when accessible; unrelated numerical matches were discarded.

### 5.1 Direct problem and mechanism searches

| Backend | Exact query | Inspection and consequence |
|---|---|---|
| Engine 1, recent coverage | Banerjee Chaudhury Schur stability counterexample | Banerjee–Chaudhury and Patel; the unrestricted claim collided. |
| Engine 1, recent coverage | Sinha kernel denoisers convergence counterexample | Sinha et al.; checked native kernel assumptions and norm-versus-spectrum distinction. |
| Engine 2 | kernel denoisers counterexample | Located the older kernel convergence/divergence literature. |
| Engine 2 | plug-and-play arbitrarily small convergence | No verified restricted Gaussian fixed-parameter collision found. |
| Engine 2 | kernel denoiser spectral radius counterexample | Prompted direct inspection of actual spectral claims. |
| Engine 2 | stochastic compression positive semidefinite stability | Generic compression/projection background; no exact restricted result located. |
| Engine 2 | Markov aggregation negative eigenvalues reversible | Aggregation terminology; no verified same-object theorem found. |
| Engine 2 | diagonally symmetrizable positive definite sum stability | Adjacent matrix-stability results; not an identity with the present iteration. |
| Engine 2 | oblique positive definite negative eigenvalues projection | The general loss-of-positivity mechanism is known background. |
| Engine 2 | Petrov Galerkin projection positive definite matrix unstable reduced system | Stability-preservation literature; does not establish the structured Gaussian classification. |

### 5.2 Rerun after the three-state quotient simplification

| Backend | Exact query | Inspection and consequence |
|---|---|---|
| Engine 2 | "stochastic" "completion" "quotient" "doubly stochastic" matrix | No verified special three-state identity located. |
| Engine 2 | "row normalization" "symmetric" "three" "contraction" matrix | No verified special identity located. |
| Engine 2 | "3x3" "stochastic" "quotient" norm | Adjacent stochastic/consensus material. |
| Engine 1 | stochastic matrix "rank one" "doubly stochastic" quotient norm | General completion background, no claimed novelty for that tool. |
| Engine 1 | symmetric row normalized "three" "Euclidean" contraction | No verified subsuming result located. |
| Engine 1 | "kernel" "three-dimensional" "Schur stability" | No verified restricted dimension result located. |
| Engine 2 | stochastic matrix second singular value centered matrix contraction coefficient | General centring and disagreement-norm aliases. |
| Engine 2 | row stochastic matrix additive rank one doubly stochastic projection | Johnson's similarity paper identified as relevant adjacent literature. |
| Engine 2 | three by three normalized kernel matrix nonexpansive Euclidean seminorm | No verified collision located. |
| Engine 2 | "Row Stochastic Matrices Similar to Doubly Stochastic Matrices" | Johnson (1981) abstract inspected; full text unavailable. |
| Engine 2 | "Euclidean" "seminorm" "three" consensus | No exact Gaussian/unit-diagonal theorem identified. |
| Engine 2 | "normalized" "kernel" "doubly stochastic" "three" contraction | No exact match located. |
| Engine 2 | "unit diagonal" "three" "stochastic" "contraction" | No verified special contraction theorem located. |
| Engine 2 | "unit diagonal" "row stochastic" "doubly stochastic" | No verified same-object result located. |
| Engine 2 | "doubly stochastic" "a(1+c)" | Pham–Spitkovsky's three-state numerical-range paper surfaced; screened as adjacent. |
| Engine 1 | "On 3-by-3 row stochastic matrices" | Publisher metadata and author-paper extract concern numerical ranges and unitary reducibility. |

### 5.3 Exact fingerprints, convolution, gap law and follow-ups

| Backend | Exact query | Inspection and consequence |
|---|---|---|
| Engine 2 | "329" "309" "109" "253" "stochastic" | Broad numerical coincidences; no verified matrix match. |
| Engine 2 | "-8 -8 17 -1" "kernel" | No verified vector/construction match. |
| Engine 2 | "1973" "74272809" "236896251" "4748139" | No verified cubic match. |
| Engine 2 | "1973" "25000000" "kernel" | Unrelated numerical/kernel hits. |
| Engine 2 | "812500" "133" "plug-and-play" | No verified threshold-constant match. |
| Engine 2 | "31/26" "133/812500" | No verified perturbation-formula match. |
| Engine 2 | "13/1272" "241/5754" | No verified compression match. |
| Engine 2 | "6287/8364672" | No verified circulant coefficient match. |
| Engine 2 | "W" "e(e" "doubly stochastic" "quotient" | Broad quotient/consensus hits, no special identity match. |
| Engine 2 | "W+e(e^T-e^TW)/3" | Formula indexing weak; no verified match. |
| Engine 2 | "W+e(eᵀ−eᵀW)/3" | Formula indexing weak; no verified match. |
| Engine 1 | "Gaussian" "circulant" "plug-and-play" "divergence" | No verified structured convolution counterexample located. |
| Engine 2 | "Gaussian" "kernel denoiser" "positive" "counterexample" | Returned the known 2021 paper; exact restrictions checked there. |
| Engine 2 | "plug-and-play" "mixing gap" "stability" | Mostly unrelated terminology; generic spectral-gap importance is already in Sinha et al. |
| Engine 2 | "Study on Convergence of Plug-and-Play ISTA With Kernel Denoisers" | Corrected the title to “With Adaptive-Kernel Denoisers.” |
| Engine 2 | "9534481" "convergence" | Publisher abstract and DOI verified. |
| Engine 1 | "Liu" "Xing" "Sun" "3111594" pdf | No accessible complete primary full text obtained. |
| Engine 1 | "Study on Convergence" "Adaptive-Kernel" filetype:pdf | References located; no complete primary text obtained. |
| Engine 2 | "adaptive-kernel denoisers" "assumption" | Conditional convergence abstract, without the exact general hypothesis. |
| Engine 2 | "Study on Convergence" "Adaptive-Kernel" "Theorem" | No additional accessible theorem text resolving that hypothesis. |

Literal coefficient searches are weak negative evidence because mathematical formulae may not be indexed accurately. Searches were supplemented with the structural aliases and direct theorem inspection.

## 6. Verified collisions and non-collisions

### 6.1 Patel: verified unrestricted collision

The posted response uses

\[
B=\frac1{11}\begin{pmatrix}9&-4\\-4&3\end{pmatrix},
\qquad W_s=\begin{pmatrix}0&1\\v&1-v\end{pmatrix},
\quad v=\frac1{5-s},\quad0<s<2.
\]

Its stationary vector is \((v,1)/(1+v)\), and
\(\pi_s^{\mathsf T}Be=s/[11(6-s)]>0\). Direct algebra gives

\[
\det(I-P_t)=-\frac{t(t-s)}{11(5-s)},\qquad
\det(I-R_t)=-\frac{t(t-s)}{(5-s)(t+1)(t+11)}.
\]

Both are negative for \(s<t<2\). Thus the posted response already has a vanishing failure threshold. The displayed matrices and determinant formulae were checked independently in this session, rather than accepted on the platform's verification label. Its signed Hessian, negative non-Perron eigenvalue of \(W_s\), and varying stationary vector distinguish it from the frozen Gaussian claim. These distinctions justify a restricted contribution, not a second claim to the broad disproof. [Patel, n.d.]

### 6.2 Gavaskar et al. 2021: actual earlier ADMM divergence

Appendix C gives a signed measurement vector \(a^{\mathsf T}=(0.8295,-0.5586)\), a positive kernel with diagonal entries \(0.1102,0.3774\), and penalty one. The source establishes an eigenvalue greater than one for the lifted standard ADMM matrix. The lifted-to-reduced identity in Section 2 translates it exactly into the present \(R_t\) convention. This is spectral divergence, not merely failure of a convenient norm.

The example does not have a unit-diagonal Gaussian affinity or a positive stochastic forward operator. Those exclusions are explicit properties of its matrices. The present work must nevertheless credit it as an older kernel ADMM counterexample. [Gavaskar et al., 2021, Appendix C.]

### 6.3 Sinha et al. 2025: native assumptions, a norm obstruction, and scaled algorithms

Assumption 1 requires a nonnegative irreducible positive semidefinite kernel with unit diagonal. Remark 1 explicitly records the unresolved standard-ISTA convergence question. Section II-C supplies a Gaussian example showing failure of contraction in the prescribed weighted norm; its construction does not show an eigenvalue outside the unit disk. The subsequent convergence theorems concern modified, scaled iterations.

The source already discusses the role of the denoiser's spectral gap in convergence rates. Accordingly, “the gap matters” is not a novelty claim here. The narrower contribution is a uniform stability bound with a family proving that linear order in the gap cannot be improved under the stated fixed-parameter restrictions. [Sinha et al., 2025.]

### 6.4 Banerjee–Chaudhury and earlier positive theorems

Banerjee–Chaudhury's Conjecture 1.1 is the broad statement already contradicted by Patel. Their Theorem 1.1 is local for each fixed instance and is compatible with the shrinking thresholds here. Their symmetric and doubly stochastic cases do not automatically apply to a row-normalised nonsymmetric Gaussian kernel. Gavaskar–Chaudhury (2020) also gives a sufficiently small, instance-dependent interval under its assumptions. [Banerjee & Chaudhury, 2025; Gavaskar & Chaudhury, 2020.]

Liu, Xing and Sun (2021) is additionally relevant. Its publisher-indexed abstract states conditional convergence for adaptive kernels and extends the inpainting step range to \((0,2)\). Full text could not be retrieved: the article page returned a bot-verification wall, and the PDF endpoint returned an access error. The exact additional hypothesis was therefore not audited. The later Sinha paper cites this work while retaining the general unresolved question, which is useful context but is not a substitute for checking Liu et al.'s full theorem.

### 6.5 Adjacent matrix results

Johnson (1981) addresses when a row stochastic matrix is similar to a doubly stochastic one. Its abstract was accessible; the full paper was not. The abstract does not establish the exact three-state positivity identity used here. Excluding a stronger hidden collision from this paper is provisional.

Pham and Spitkovsky (2023) study numerical-range shapes, Kippenhahn curves and unitary reducibility of three-state row stochastic matrices. The abstract and accessible author-paper extract did not state the unit-diagonal row-normalisation contraction result. Publisher full-text retrieval was inconsistent and ultimately failed, so this is a targeted screening, not a complete exclusion.

The general fact that oblique reduction can destroy positivity or stability is established background in projection/model-reduction literature. No originality is assigned to that mechanism alone.

## 7. Missingness and coverage audit

The search used live web retrieval, two search backends, arXiv full texts, publisher records, author publication lists, MathDB, and selected matrix-analysis/model-reduction aliases. Direct full texts were inspected for Banerjee–Chaudhury, Sinha et al., Gavaskar–Chaudhury, Gavaskar et al., and Patel's response. Recent author-page and arXiv follow-ups were also screened; no verified exact restricted collision was identified.

Material limitations remain:

- Johnson (1981) and Liu et al. (2021) were not available in full text through the available route. Pham–Spitkovsky's publisher retrieval was incomplete.
- Formula search engines may discard signs, denominators, transposes, or mathematical formatting. Exact-number searches have low recall.
- No exhaustive MathSciNet, zbMATH, Scopus or Web of Science search was performed, and no claim of comprehensive citation chaining is made.
- Unindexed papers, non-English literature, unpublished work, theses and very recent postings can be missed.
- No author was contacted. No external specialist has certified novelty.
- No OEIS audit is claimed: the object under study is a finite family of matrices and determinant polynomials, not a natural counting sequence.
- The exact verifiers validate arithmetic and witnesses, not publication priority.

These limitations justify BOUNDED UNCERTAINTY rather than an unqualified “new theorem” or “first counterexample” declaration.

## 8. Outcome by unit of contribution

| Unit | Outcome | Permissible claim |
|---|---|---|
| General disproof of the Banerjee–Chaudhury interval | COLLISION | Cite Patel. |
| Arbitrarily small failure thresholds without the Gaussian/fixed-parameter restrictions | COLLISION | Cite Patel; do not count this as new. |
| General standard kernel ADMM divergence | COLLISION | Cite Gavaskar et al. 2021. |
| Native Gaussian four-state family with fixed positive \(A,\pi\) and uniform \(K\)-conditioning | BOUNDED UNCERTAINTY | Explicitly proved here; no subsuming result located in the inspected corpus. |
| Full stability in dimensions at most three and sharp dimension four | BOUNDED UNCERTAINTY | Candidate structural theorem; adjacent matrix full-text check remains desirable. |
| Gap-dependent sufficient bound plus matching-order obstruction | BOUNDED UNCERTAINTY | Claim sharp order within the stated normalisation; do not claim a sharp constant or the first recognition of gap dependence. |
| Positive circulant Gaussian family | BOUNDED UNCERTAINTY | Candidate extension to positive periodic convolution; no minimal circulant dimension claimed. |
| Determinant identity, quotient norms, Gaussian distance calculations, oblique compression | STANDARD TOOLS / BRIDGE | Explanatory tools and an explicit bridge between literatures, not separate discoveries. |
| Exact finite and polynomial certificates | VERIFIED COMPUTATION | Reproducibility evidence; this status is separate from novelty. |

## 9. Claim wording and next decision

The manuscript's present framing is appropriate:

> We prove a sharp dimension threshold within the constant-preserving Hessian class and construct genuine Gaussian-kernel examples with fixed positive forward model, fixed stationary weights and uniformly bounded affinity conditioning. A separate example uses a positive circulant forward model. A bounded search did not locate these restricted results; broader counterexamples are already known.

This is a paraphrase of the permitted contribution claim, not a quotation from an external source.

The mathematics is sufficiently concrete to circulate for specialist assessment. Publication priority should remain provisional. The most valuable remaining literature work is obtaining Johnson's and Liu et al.'s complete papers and asking a matrix-analysis/inverse-problems specialist to check whether the three-state identity or the structured Gaussian obstruction is already implicit in a stronger theorem. These are confidence improvements, not presently identified defects in the proofs.

## 10. Rerun triggers

Rerun this gate if the claim is broadened beyond \(Be=b e\), the kernel is reinterpreted as image patches, a conventional Gaussian convolution is claimed, the gap constant is declared optimal, a smaller circulant construction is found, or a new source is reported to contain a related quotient-contraction or Gaussian divergence result. A mere change of notation, vertex order, or Hessian units is not a new object.

## References

Banerjee, D., & Chaudhury, K. N. (2025). *On the Schur stability of some image reconstruction operators*. arXiv. https://arxiv.org/html/2509.06715v1

Gavaskar, R. G., & Chaudhury, K. N. (2020). Plug-and-play ISTA converges with kernel denoisers. *IEEE Signal Processing Letters, 27*, 610–614. https://doi.org/10.1109/LSP.2020.2986643. Full text: https://arxiv.org/html/2004.03145v2

Gavaskar, R. G., Athalye, C. D., & Chaudhury, K. N. (2021). On plug-and-play regularization using linear denoisers. *IEEE Transactions on Image Processing, 30*, 4802–4813. https://doi.org/10.1109/TIP.2021.3075092. Full text: https://arxiv.org/html/2105.05177v1

Johnson, C. R. (1981). Row stochastic matrices similar to doubly stochastic matrices. *Linear and Multilinear Algebra, 10*(2), 113–130. https://doi.org/10.1080/03081088108817402 [Abstract inspected; full text unavailable.]

Liu, T., Xing, L., & Sun, Z. (2021). Study on convergence of plug-and-play ISTA with adaptive-kernel denoisers. *IEEE Signal Processing Letters, 28*, 1918–1922. https://doi.org/10.1109/LSP.2021.3111594. Publisher record: https://ieeexplore.ieee.org/document/9534481/ [Publisher-indexed abstract inspected; full text unavailable.]

Patel, S. (n.d.). *Counterexample response to “Schur stability conjecture for image reconstruction operators”* [Online mathematical response]. MathDB. Retrieved 20 October 2026, from https://mathdb.com/p/370245/schur-stability-conjecture-for-image-reconstruction-operator

Pham, N., & Spitkovsky, I. M. (2023). On 3-by-3 row stochastic matrices. *Special Matrices, 11*(1), Article 20230103. https://doi.org/10.1515/spma-2023-0103 [Abstract and accessible author-paper extract screened; complete publisher retrieval unsuccessful.]

Sinha, A., Kumar, B., Athalye, C. D., & Chaudhury, K. N. (2025). Linear convergence of plug-and-play algorithms with kernel denoisers. *IEEE Transactions on Signal Processing, 73*, 2646–2659. https://doi.org/10.1109/TSP.2025.3573044. Full text: https://arxiv.org/html/2505.15318v1

