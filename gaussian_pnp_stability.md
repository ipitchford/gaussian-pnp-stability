---
title: "Gaussian-kernel plug-and-play iterations can be unstable"
subtitle: "A sharp dimension threshold for irreducible unit-diagonal kernels, the role of the mixing gap, and why unit diagonals matter"
author: "Anonymous"
date: "27 October 2026 · Version 1.0.0-candidate (unrefereed) · DOI [10.5281/zenodo.23690125](https://doi.org/10.5281/zenodo.23690125) · Code and data: <https://github.com/ipitchford/gaussian-pnp-stability>"
lang: en-GB
---

## Abstract

Plug-and-play (PnP) methods for linear inverse problems alternate a data step with a denoiser. When the denoiser is a row-normalised kernel smoother, such as a Gaussian non-local-means filter, it is not known in general whether the iteration converges for every step size that is admissible for the data term. Sinha et al. (2025) state that they had found neither a proof nor a counterexample. We answer the question for the two standard linear iterations (PnP-ISTA and PnP-ADMM).

The answer depends on the dimension. With irreducible unit-diagonal kernels and a data Hessian that preserves constants, both iterations are stable in dimensions one to three, for every admissible step. In dimension four they need not be: we give explicit Gaussian kernels and a fixed positive forward operator for which both iterations diverge at every step up to one half. A twelve-dimensional example shows the same failure for a positive convolution (deblurring) operator, but only for very small steps. If unit diagonals are not required, failure already occurs in dimension three.

The cause is a denoiser whose similarity graph is almost disconnected and averages unequally. The divergence is extremely slow, about $1.7\times10^{-5}$ per iteration at most, and the affinities involved are tiny. The practical lesson is that safe step sizes must scale with the denoiser's mixing gap, or the denoiser should be symmetrised. We prove both remedies.

## 1. Introduction

### 1.1 The question

A linear inverse problem seeks an unknown vector $x$ from observations of $Ax$. For a least-squares data term the Hessian is $B=A^{\mathsf T}A$. A plug-and-play method alternates a step on this data term with a denoiser. For a fixed linear denoiser $W$, the two standard iteration matrices are

\[
P(t)=W(I-tB),\qquad
R(t)=I-W+(I+tB)^{-1}(2W-I),\qquad t>0.
\tag{1}
\]

$P(t)$ is the iteration matrix of PnP-ISTA (proximal-gradient) and $R(t)$ that of PnP-ADMM (Banerjee and Chaudhury, 2025, Section 2). The iterations converge for every starting point exactly when these matrices are *Schur stable*, that is, when every eigenvalue has modulus less than one. A real eigenvalue greater than one gives a direction in which errors grow geometrically.

For gradient descent without a denoiser, every step $0<t<2/\rho(B)$ is safe. Kernel denoisers, the classical non-local-means and bilateral filters, are linear once their weights are fixed:

\[
W=D^{-1}K,\qquad D=\operatorname{diag}(Ke),
\tag{2}
\]

where $K$ is a symmetric nonnegative affinity matrix and $e$ is the all-ones vector. Gaussian affinities have the form

\[
K_{ij}=\exp(-\|x_i-x_j\|_2^2)
\tag{3}
\]

for feature vectors $x_i$. Row normalisation makes $W$ an averaging operator but destroys symmetry. Whether the interval $0<t<2/\rho(B)$ remains safe for such $W$ is open. Banerjee and Chaudhury (2025, Conjecture 1.1) conjectured that it is, under conditions on $B$. Sinha et al. (2025; arXiv version, Section II-A, in the paragraph following Remark 1) observe that PnP-ISTA with kernel denoisers seems to converge in practice, but that they "have not yet found a formal proof or a counterexample". This paper settles that question for unit-diagonal kernels.

### 1.2 The answer in brief

Throughout, $B$ is symmetric positive semidefinite and preserves constants, $Be=be$ with $b>0$, and the kernel is *irreducible* (its similarity graph is connected) with unit diagonal and entries at most one.

- **Dimensions one to three are always stable** (Theorem 1(a)), for every $0<t<2/\rho(B)$. This reduces to Banerjee and Chaudhury's Theorem 1.3 by a new positivity identity.
- **Dimension four can be unstable** (Theorem 1(b)). An explicit Gaussian kernel makes both iterations diverge at any prescribed step $0<t\leq1/2$, with a fixed positive forward operator. The exact boundary of this family is $t^*\approx0.5313$ (Theorem 4).
- **The cause is a small mixing gap.** Every iteration is stable for $t\lesssim g/\sqrt{\kappa_\pi}$ (Theorem 2), and in our family the first instability occurs at about $6109\,g$, so linear scaling in $g$ cannot be improved.
- **Convolution operators are not immune** (Theorem 3), but the certified failure window is only $0<t<1.9\times10^{-4}$.
- **Unit diagonals matter** (Theorem 5). Without them, failure occurs already in dimension three.
- **Symmetric denoisers are stable in every dimension**, exactly when no nonzero vector is fixed by $W$ and annihilated by $B$ (Theorem 6). Irreducible symmetrically normalised kernels always meet this condition when $Be=be$ (Corollary 1). The spectral-norm argument is due to Athalye, Chaudhury and Kumar (2023).

### 1.3 Vocabulary

- A **stochastic** matrix has nonnegative entries and row sums one.
- A nonnegative matrix is **irreducible** if its directed support graph is strongly connected, and **primitive** if some power of it is strictly positive. For the symmetric kernels used here, and for reversible matrices with symmetric support, irreducibility is the same as connectedness of the undirected similarity graph.
- A stochastic $W$ has a **stationary distribution** $\pi$, a probability vector with $\pi^{\mathsf T}W=\pi^{\mathsf T}$. It is **reversible** if $\pi_iW_{ij}=\pi_jW_{ji}$; every $W$ of the form (2) is reversible, with $\pi\propto Ke$. We write $\kappa_\pi=\max_i\pi_i/\min_i\pi_i$ for the imbalance of the averaging weights.
- The **mixing gap** of a reversible $W$ with real eigenvalues $1=\lambda_1\geq\lambda_2\geq\cdots$ is $g=1-\lambda_2$. It measures how quickly averaging spreads across the similarity graph. It is small when the graph is almost disconnected, and it is zero when the graph is disconnected.
- **The quotient modulo constants.** Constants form an invariant subspace of $W$ and of $B$. The quotient keeps only the differences between coordinates; all instabilities live there.
- An **oblique compression** of $B$ uses one map to send group levels to vertices and a different, weighted, map to average them back. Unlike an orthogonal compression, it can turn a positive definite $B$ into a matrix with a negative eigenvalue.
- The **spectral condition number** of a symmetric positive definite $M$ is $\kappa_2(M)=\lambda_{\max}(M)/\lambda_{\min}(M)$.

### 1.4 Summary of results

| Kernel class | Forward / Hessian class | Smallest unstable dimension | Unstable steps certified | Source |
|---|---|---|---|---|
| any primitive $W$; signed $B$ with $Be\neq\rho(B)e$ | signed | 2 | thresholds $\to0$ | Patel (online post, see §2) |
| non-unit diagonal, signed forward vector (ADMM only) | signed | 2 | some $t$ | Gavaskar et al. (2021, App. C) |
| irreducible, unit diagonal, $K_{ij}\le1$ | $B\succeq0$, $Be=be$ | **4** (none for $n\le3$) | $0<t<t^*\approx0.5313$, varying kernel | Thms 1, 4 |
| same, Gaussian, positive circulant $A$ | convolution | $\le12$ | $0<t<1.9\times10^{-4}$ on $\epsilon=t^2$ | Thm 3 |
| reversible positive kernel, unequal diagonal | positive, fixed | **3** (none for $n\le2$) | an interval $(t_-,t_+)$ | Thm 5, Prop. 1 |
| symmetric $W$, spectrum in $(-1,1]$, $\ker(I-W)\cap\ker B=\{0\}$ | $B\succeq0$ | never unstable | none | Thm 6; Athalye et al. (2023) |

### 1.5 Roadmap

Section 2 places the work in the literature. Section 3 states the main theorem. Section 4 explains the mechanism with almost no algebra. Section 5 proves stability in dimensions one to three. Sections 6 and 7 build and certify the four-dimensional counterexample. Section 8 gives the step-size guarantee in terms of the mixing gap, and a complete two-parameter classification of the four-state family that settles where stability is first lost. Section 9 treats convolution operators. Section 10 shows that the unit-diagonal assumption is essential. Section 11 characterises when symmetric denoisers are stable, extending Athalye et al. (2023). Section 12 explains what the results mean in practice. Section 13 collects scope limits and non-claims, and Section 14 lists open problems. Reproduction details, the large constants and the certificates for Theorems 3 and 4 (Appendix C) are at the end.

Figure 1 summarises the four-state family and the three-state family at a glance.

![(a) The exact unstable region of the four-state family in the plane of weak-link strength $\epsilon$ and step $t$ (shaded, Theorem 4). The dashed curve is the one-parameter family of Theorem 1(b); it leaves the region at $t^*\approx0.5313$. For $\epsilon>\epsilon_{\rm crit}\approx1.69\times10^{-5}$ there is no instability; at $\epsilon=\epsilon_{\rm crit}$ there is a single neutral step with eigenvalue one. (b) Excess spectral radius per iteration along that family, for PnP-ISTA and PnP-ADMM. The curves coincide to plotting accuracy, and the excess never exceeds $1.72\times10^{-5}$. (c) The true first unstable step against the mixing gap, compared with the guaranteed safe step of Theorem 2; both are linear in $g$, about 8,600 apart. (d) The exact unstable region of the three-state family with unequal diagonal entries (Theorem 5): stability is lost and then regained.](gaussian_stability_diagram.pdf){width=100%}

## 2. Prior work and positioning

**The conjecture and its refutation.** Banerjee and Chaudhury (2025) prove that each fixed primitive stochastic $W$ is stable on some sufficiently small interval when $\pi^{\mathsf T}Be>0$ (their Theorem 1.1), and that the full interval $0<t<2/\rho(B)$ is safe in several special cases. These are symmetric $W$ (Theorem 1.2, which they adapt from Athalye, Chaudhury and Kumar, 2023), doubly stochastic $W$ with $W^{\mathsf T}W$ and $WW^{\mathsf T}$ irreducible (Theorem 1.3), diagonal $B$ (Proposition 1.1) and $B=\alpha I+\beta ee^{\mathsf T}$ (Proposition 1.2). Their Conjecture 1.1 extends this to any primitive stochastic $W$ and positive semidefinite $B$ with $Be\leq\rho(B)e$ and $\pi^{\mathsf T}Be>0$.

A counterexample to that conjecture was posted by S. Patel on the MathDB problem site (https://mathdb.com/p/370245/schur-stability-conjecture-for-image-reconstruction-operator, accessed 20 October 2026). The post has no Crossref or arXiv record, so it is acknowledged here rather than listed among the references. MathDB labels it an unverified claim; its determinant identities are elementary and were checked by a referee of this paper. Patel's construction is two-dimensional. Its denoiser $W=\begin{pmatrix}0&1\\ v&1-v\end{pmatrix}$ has a zero diagonal entry, so it cannot come from a kernel with positive diagonal. Its Hessian $B=\frac1{11}\begin{pmatrix}9&-4\\-4&3\end{pmatrix}$ has negative entries and does not preserve constants. In that family the unstable thresholds tend to zero. We do not claim the first disproof of Conjecture 1.1.

**Kernel denoisers specifically.** Gavaskar, Athalye and Chaudhury (2021, Appendix C) showed that standard PnP-ADMM can diverge with a two-dimensional kernel denoiser; their example uses a signed forward vector and unequal diagonal entries. Gavaskar and Chaudhury (2020) proved PnP-ISTA convergence with kernel denoisers for sufficiently small, instance-dependent steps. Liu, Xing and Sun (2021) report conditional convergence for adaptive-kernel denoisers; we were able to read only its abstract, not the full paper, so its exact hypotheses have not been checked against ours. Sinha et al. (2025) prove linear convergence for *scaled* PnP updates, show that kernel PnP-ISTA need not contract in a particular norm (which does not imply spectral instability), and pose the question answered here. Our matrices are the unscaled ones in (1). Our kernels satisfy their standing assumptions: the four-state kernels are symmetric, unit-diagonal and positive definite.

**The main convergence routes for PnP.** Most convergence theory for PnP assumes something the row-normalised kernel denoiser lacks:

- nonexpansive, averaged or firmly nonexpansive denoisers (Ryu et al., 2019; Pesquet et al., 2021);
- denoisers that are gradient steps or proximal maps of a potential (Hurault, Leclaire and Papadakis, 2021, 2022);
- a symmetric Jacobian, as in regularisation by denoising (Romano, Elad and Milanfar, 2017);
- bounded denoisers with continuation (Chan, Wang and Elgendy, 2017);
- symmetric, doubly stochastic smoothing matrices (Sreehari et al., 2016).

A row-normalised kernel denoiser is not symmetric, need not be nonexpansive in the Euclidean norm, and is not a proximal map. It therefore falls outside all of these routes. This is why the question is genuinely open for kernel denoisers, and why a counterexample within this class is informative.

**What is new here.** As far as we know, the following are new: the sharp dimension four for irreducible unit-diagonal kernels when $B$ preserves constants; the explicit Gaussian family with a fixed positive forward operator and its exact classification; the order-$g$ obstruction; the certified positive-circulant family; the three-state reentrant classification; and, for symmetric denoisers, the extension of Athalye et al.'s (2023) characterisation to Schur stability of denoisers that need not be stochastic or positive semidefinite, with Corollary 1 for the symmetric normalisation of irreducible kernels. The low-dimensional positive result is, at heart, Banerjee and Chaudhury's Theorem 1.3. Our contribution there is the positivity identity (6)–(7) that brings three-state kernels within its scope. Section 13 records the limits of our literature search.

## 3. Main theorem

**Theorem 1 (sharp dimension for irreducible unit-diagonal kernels).**

**(a) Stability in dimensions one to three.** Let $n\leq3$. Suppose $K\in\mathbb R^{n\times n}$ is symmetric, nonnegative and irreducible, with $K_{ii}=1$ and $K_{ij}\leq1$, and let $W$ be given by (2). Suppose $B$ is symmetric positive semidefinite with $Be=be$ and $b>0$. Then both matrices in (1) are Schur stable for every $0<t<2/\rho(B)$. Positive semidefiniteness of $K$ is not required. For a nonnegative unit-diagonal positive semidefinite kernel, the entry bound $K_{ij}\le1$ holds automatically.

**(b) Failure in dimension four.** There are a fixed strictly positive symmetric stochastic positive definite matrix $A\in\mathbb R^{4\times4}$ and a fixed positive probability vector $\pi$ such that, for every prescribed $0<t\leq1/2$, there is a four-point Gaussian kernel $K_t$ of the form (3) with the following properties:

- $B=A^{\mathsf T}A$ is fixed, strictly positive, positive definite and stochastic, with $\rho(B)=1$ and $\kappa_2(B)=2500$;
- $W_t=\operatorname{diag}(K_te)^{-1}K_t$ is strictly positive and reversible with stationary distribution $\pi$;
- $K_t$ is positive definite with $\kappa_2(K_t)<19.22$, and all eigenvalues of $W_t$ are positive;
- both $P_t=W_t(I-tB)$ and $R_t=I-W_t+(I+tB)^{-1}(2W_t-I)$ have a real eigenvalue strictly greater than one.

**Consequently, among irreducible unit-diagonal kernels with entries at most one, and Hessians with $Be=be$, four is the smallest dimension in which instability can occur.** In particular four is minimal in the strictly positive Gaussian subclass.

Four remarks fix the scope.

*Irreducibility is essential.* Without it the threshold fails already in dimension two. Take $K=W=I_2$ and $B=ee^{\mathsf T}/2$, which is symmetric, positive semidefinite and stochastic with $\rho(B)=1$. For every $0<t<2$, $P(t)=I-tB$ has eigenvalues $1-t$ and $1$, and $R(t)=(I+tB)^{-1}$ has eigenvalues $(1+t)^{-1}$ and $1$. Neither is Schur stable.

*The data Hessian must preserve constants.* Part (a) assumes $Be=be$. Conjecture 1.1 allows the weaker condition $Be\leq\rho(B)e$, and we do not know whether part (a) survives under it (Section 14).

*The kernel depends on the step.* For each prescribed step we select a kernel and then hold it fixed during the iteration. We do not claim that one fixed kernel is unstable at all small steps; by Banerjee and Chaudhury's Theorem 1.1, it cannot be.

*The certified range is $t\le1/2$, and the exact range is $t<t^*$.* The family of part (b) is unstable exactly for $0<t<t^*\approx0.5313$ and is stable again for $t^*<t<2$ (Theorem 4). Whether some other four-point unit-diagonal Gaussian kernel is unstable at a step in $(t^*,2)$ is open (Section 14). The minimality statement is therefore matched on $(0,1/2]$ only.

## 4. The mechanism: unequal averaging on an almost disconnected graph

Imagine a denoiser whose similarity graph has two groups joined by a very weak link of strength $\epsilon$. Vectors that are nearly constant on each group change slowly under denoising. What matters is what the data step does to the *difference* between the two group levels.

A symmetric positive definite $B$ damps every direction in Euclidean geometry. The denoiser, however, averages each group with its own stationary weights. When those weights are unequal, the effective two-group data operator is an oblique compression: the map that averages a group is not the transpose of the map that assigns it a constant value. Such a compression need not be positive definite. In our example its nontrivial eigenvalue is $c=-133/812500$ (Section 8.4), which is negative.

Two effects then compete on the slow difference mode. The weak link damps it by an amount proportional to the mixing gap, about $(31/26)\epsilon$. The negative effective data coefficient amplifies it by about $|c|t$. Instability results when amplification wins.

**The example is a controlled perturbation of the disconnected case.** At $\epsilon=0$ the four-state kernel $K_0$ is reducible: the fourth point is isolated. Its denoiser $W_0$ has no damping of the difference mode at all, and $P$ and $R$ are already unstable for every $0<t<1330000/2254167\approx0.590$ (Theorem 4 at $\epsilon=0$). Theorem 1(b) attaches the isolated point by affinities of order $\epsilon=t/100000$. This makes the kernel strictly positive and Gaussian while keeping the reducible instability, with exact uniform control. The novelty lies in that control (a fixed forward operator, a fixed stationary distribution, bounded conditioning, genuine Gaussian geometry and an exact certificate), not in a new limiting mechanism. Practical kernels that truncate small affinities sit at, not near, the reducible limit (Section 12).

Why four? On two vertices a connected unit-diagonal kernel has equal row sums, so its averaging is uniform. A non-uniform averaging block first becomes possible on three vertices, and a second group needs a fourth. This is only intuition. The proof of Theorem 1(a) is stronger: it establishes stability for every admissible kernel in dimensions up to three, with no assumption on the graph structure.

## 5. Why dimensions one, two and three are stable

We first separate constants from differences. Let $U$ have orthonormal columns spanning $e^\perp$, and put

\[
Q=U^{\mathsf T}WU,\qquad C=U^{\mathsf T}BU.
\]

In the orthonormal basis made of $e/\sqrt n$ and the columns of $U$,

\[
W\sim\begin{pmatrix}1&r^{\mathsf T}\\0&Q\end{pmatrix},
\qquad
B\sim\begin{pmatrix}b&0\\0&C\end{pmatrix}.
\tag{4}
\]

$Q$ describes the action of $W$ on vectors modulo constants. Although $W$ need not preserve $e^\perp$, (4) represents this quotient action correctly. The constant-mode eigenvalues of $P$ and $R$ are $1-tb$ and $1/(1+tb)$. Since $b$ is an eigenvalue of $B$, $b\le\rho(B)$, and both are stable for $0<t<2/\rho(B)$. It remains to control the quotients. All norm arguments below are applied to the complexified quotient space, so they cover complex eigenvalues.

### 5.1 An elementary three-state positivity identity

Write

\[
K=\begin{pmatrix}1&k_{12}&k_{13}\\k_{12}&1&k_{23}\\k_{13}&k_{23}&1\end{pmatrix},
\quad
d_1=1+k_{12}+k_{13},\quad d_2=1+k_{12}+k_{23},\quad d_3=1+k_{13}+k_{23},
\]

where $0\leq k_{12},k_{13},k_{23}\leq1$. Define

\[
S=W+\frac13e(e^{\mathsf T}-e^{\mathsf T}W).
\tag{5}
\]

This rank-one correction gives $S$ row and column sums one, and it does not change the quotient: $U^{\mathsf T}SU=Q$, because $U^{\mathsf T}e=0$. Its essential property is positivity. For example,

\[
\begin{aligned}
3S_{21}
&=1-\frac1{d_1}+\frac{2k_{12}}{d_2}-\frac{k_{13}}{d_3}\\
&=\frac{k_{12}(1+k_{23})+k_{13}k_{23}}{d_1d_3}+\frac{2k_{12}}{d_2}\geq0.
\end{aligned}
\tag{6}
\]

For a diagonal entry,

\[
3S_{11}=1+\frac2{d_1}-\frac{k_{12}}{d_2}-\frac{k_{13}}{d_3}
\geq\frac2{d_1}>0,
\tag{7}
\]

because

\[
\frac{k_{12}}{d_2}+\frac{k_{13}}{d_3}
\leq\frac{k_{12}}{1+k_{12}}+\frac{k_{13}}{1+k_{13}}\leq1.
\]

The last inequality is equivalent to $k_{12}k_{13}\leq1$. This is exactly where the hypothesis $K_{ij}\le1$ is used. The general entry $S_{ij}$ is not symmetric in $(i,j)$. The remaining entries follow from (6) and (7) by applying every permutation of the three vertices, which permutes the $k$'s and covers both ordered roles $(i,j)$ and $(j,i)$. Connectivity makes every off-diagonal entry strictly positive: either the direct edge is positive, or the other two edges give a positive two-edge term in (6).

### 5.2 Reduction to Banerjee and Chaudhury's Theorem 1.3

Thus $S$ is strictly positive and doubly stochastic, so $S^{\mathsf T}S$ and $SS^{\mathsf T}$ are strictly positive and hence irreducible. Banerjee and Chaudhury's Theorem 1.3 therefore applies to $S$ and $B$ (note $Be=be\neq0$), so $S(I-tB)$ and $I-S+(I+tB)^{-1}(2S-I)$ are stable for $0<t<2/\rho(B)$. Their spectra coincide with those of $P(t)$ and $R(t)$, because $W$ and $S$ have the same quotient $Q$ and the same constant mode. This proves Theorem 1(a) for $n=3$.

For completeness we include the short self-contained argument, which follows theirs. Jensen's inequality gives $\|Sx\|_2^2\leq\sum_{i,j}S_{ij}|x_j|^2=\|x\|_2^2$, with equality only for constant $x$, since every row of $S$ gives positive weight to every coordinate. By compactness of the unit sphere in $e^\perp$,

\[
\|Q\|_2<1.
\tag{8}
\]

For $0<t<2/\rho(B)$, $\|I-tC\|_2\le1$, so $\|Q(I-tC)\|_2<1$ and the quotient of $P$ is stable. For $R$ the quotient is $R_q=I-Q+(I+tC)^{-1}(2Q-I)$. Suppose $R_qx=\lambda x$ with $x\neq0$ and $|\lambda|\geq1$. Multiplying by $I+tC$ gives

\[
[\lambda I+(\lambda-1)tC]x=(I-tC)Qx.
\tag{9}
\]

For every eigenvalue $\beta\in[0,\rho(B)]$ of $C$, $|\lambda(1+t\beta)-t\beta|\geq|\lambda|(1+t\beta)-t\beta\geq1$. The left-hand matrix in (9) is therefore invertible and $x=\varphi_\lambda(C)Qx$, with $\varphi_\lambda(\beta)=(1-t\beta)/(\lambda(1+t\beta)-t\beta)$. Since $|1-t\beta|\leq1$, $\|\varphi_\lambda(C)\|_2\leq1$, and $\|x\|_2\leq\|Qx\|_2<\|x\|_2$, a contradiction.

For $n=2$, a row-normalised symmetric unit-diagonal kernel has equal row sums, so it is itself symmetric and stochastic, and the same argument applies. The case $n=1$ is scalar. This proves Theorem 1(a).

## 6. The explicit four-dimensional family

### 6.1 Fix the forward operator once

Let $J=ee^{\mathsf T}$ and $v=(-8,-8,17,-1)^{\mathsf T}$. Set

\[
A=\frac{20I+245J+vv^{\mathsf T}}{1000}
=\frac1{1000}\begin{pmatrix}
329&309&109&253\\
309&329&109&253\\
109&109&554&228\\
253&253&228&266
\end{pmatrix},\qquad B=A^2.
\tag{10}
\]

Because $e^{\mathsf T}v=0$ and $\|v\|_2^2=418$, the orthogonal spaces spanned by $e$, by $v$, and by their common orthogonal complement give

\[
\sigma(A)=\left\{1,\frac{219}{500},\frac1{50},\frac1{50}\right\},
\qquad
\sigma(B)=\left\{1,\frac{47961}{250000},\frac1{2500},\frac1{2500}\right\}.
\tag{11}
\]

The entries of $A$ are at least $109/1000$. $A$ is symmetric, stochastic and positive definite, so $B=A^{\mathsf T}A$ is strictly positive, symmetric, stochastic and positive definite, with $\rho(B)=1$ and $\kappa_2(B)=2500$. For later use,

\[
B=\frac1{2500}I+\frac{2499}{10000}J+\frac{229}{500000}vv^{\mathsf T}.
\tag{12}
\]

### 6.2 Vary only the affinity matrix

Define

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
\end{pmatrix},
\tag{13}
\]

and, for $\epsilon\geq0$,

\[
K_\epsilon=K_0+\epsilon H,\qquad
D_\epsilon=\operatorname{diag}(K_\epsilon e),\qquad W_\epsilon=D_\epsilon^{-1}K_\epsilon.
\tag{14}
\]

Theorem 1(b) uses $\epsilon=t/100000$, and we then write $K_t,W_t,P_t,R_t$. The equal row-sum vectors $K_0e=He=(2,2,6/5,1)^{\mathsf T}$ give

\[
D_\epsilon=(1+\epsilon)\operatorname{diag}(2,2,6/5,1),\qquad
\pi=\frac1{31}(10,10,6,5)^{\mathsf T}.
\tag{15}
\]

The stationary distribution is therefore exactly fixed, and $\kappa_\pi=2$. Symmetry of $K_\epsilon$ gives detailed balance. For $\epsilon>0$ every entry of $K_\epsilon$ is positive and its diagonal is one.

The eigenvalues of $K_0$ are $1/10$ (on $(1,-1,0,0)^{\mathsf T}$), $1$ (on the fourth coordinate) and $(29\pm\sqrt{89})/20$ (on vectors with $x_1=x_2$, $x_4=0$). Hence $\kappa_2(K_0)=(29+\sqrt{89})/2\approx19.217$. Since $\|H\|_2\leq2$ (by symmetry and its maximum absolute row sum), Weyl's inequality gives, for $0<\epsilon\leq1/200000$,

\[
\lambda_{\min}(K_\epsilon)\geq\frac1{10}-\frac1{100000},
\qquad \lambda_{\max}(K_\epsilon)\leq\frac{29+\sqrt{89}}{20}+\frac1{100000},
\qquad \kappa_2(K_\epsilon)<19.22.
\tag{16}
\]

The matrix $W_\epsilon$ is similar to the symmetric positive definite matrix $D_\epsilon^{-1/2}K_\epsilon D_\epsilon^{-1/2}$, so all its eigenvalues are positive. Its eigenvalue one is simple because $W_\epsilon$ is strictly positive.

## 7. An exact certificate of instability

**Strategy.** A real $4\times4$ matrix $M$ has a real eigenvalue greater than one whenever $\det(I-M)<0$: its characteristic polynomial is negative at $\lambda=1$ and positive for large $\lambda$. We compute $\det(I-P_t)$ exactly, after splitting off one symmetric direction.

Define

\[
M_t=D_t-K_t+tK_tB=D_t(I-P_t).
\tag{17}
\]

All matrices in (17) respect the exchange of coordinates 1 and 2. On the antisymmetric direction $(1,-1,0,0)^{\mathsf T}$, $M_t$ has the eigenvalue

\[
m_a=2(1+\epsilon)-\left(\frac1{10}-\frac32\epsilon\right)
\left(1-\frac{t}{2500}\right)>0,
\tag{18}
\]

since both bracketed factors lie in $(0,1)$ while the first term exceeds two. On the space with $x_1=x_2$, coordinates $(x_1,x_3,x_4)$ give a $3\times3$ matrix $\bar M_t$ (Appendix A). Exact expansion gives

\[
\det\bar M_t=t^2\left(1+\frac{t}{100000}\right)f(t),
\tag{19}
\]

where $f$ is the cubic

\[
f(t)\approx-7.89\times10^{-5}+1.49\times10^{-4}\,t-3.79\times10^{-10}\,t^2-3.04\times10^{-15}\,t^3,
\tag{20}
\]

with exact rational coefficients listed in Appendix B. The sign argument is short. The linear coefficient is less than $3/20000$ and both higher coefficients are negative. Therefore, throughout $0<t\leq1/2$,

\[
f(t)<-\frac{1973}{25000000}+\frac12\cdot\frac3{20000}
=-\frac{49}{12500000}<0.
\tag{21}
\]

So $\det M_t=m_a\det\bar M_t<0$ and, since $\det D_t>0$, $\det(I-P_t)<0$. Hence $P_t$ has a real eigenvalue greater than one.

The exact boundary of the family is the unique root $t^*\in(0,2]$ of $f$:

\[
t^*=0.53128532534\ldots,
\tag{22}
\]

so $\det(I-P_t)<0$ precisely for $0<t<t^*$, and $\det(I-P_t)>0$ for $t^*<t\le2$. Theorem 4 shows that the family is in fact Schur stable for $t^*<t<2$.

The second algorithm follows from the same certificate. Direct multiplication gives $(I+tB)(I-R_t)=I-(I-tB)W_t$, and Sylvester's determinant identity reverses the two factors:

\[
\det(I-R_t)=\frac{\det(I-W_t(I-tB))}{\det(I+tB)}
=\frac{\det(I-P_t)}{\det(I+tB)}<0.
\tag{23}
\]

The denominator is positive because $B$ is positive definite. This is a determinant identity; we do not assert that the two numerator matrices are equal.

### 7.1 The matrices are genuine Gaussian kernels

A positive definite unit-diagonal matrix need not have a Gaussian distance representation, so we verify it. Write

\[
a=\frac9{10}+\frac32\epsilon,\qquad
b_1=\frac1{10}+\frac3{10}\epsilon,\qquad
u=-\log a,\qquad w=-\log b_1,
\]

\[
h=\sqrt{w-u/4},\qquad
y=\frac{w-u/2+\log3}{2h},\qquad
z=\sqrt{-\log(\epsilon/5)-u/4-y^2},
\tag{24}
\]

and choose four feature vectors in $\mathbb R^3$:

\[
x_1=(-\sqrt u/2,0,0),\quad
x_2=(\sqrt u/2,0,0),\quad
x_3=(0,h,0),\quad
x_4=(0,y,z).
\tag{25}
\]

All square roots are real for $0<\epsilon\leq1/200000$, because $0<u<1/8$, $2<w<3$ and $\log3<2$. Hence $h^2>63/32$ and $0<w-u/2+\log3<5$, so $u/4+y^2<1/32+200/63<4<\log(10^6)\leq-\log(\epsilon/5)$. The six squared distances are

\[
\begin{aligned}
\|x_1-x_2\|_2^2&=u,\qquad
\|x_1-x_3\|_2^2=\|x_2-x_3\|_2^2=w,\\
\|x_1-x_4\|_2^2=\|x_2-x_4\|_2^2&=-\log(\epsilon/5),\qquad
\|x_3-x_4\|_2^2=-\log(3\epsilon/5).
\end{aligned}
\tag{26}
\]

For the last, substitute $h^2=w-u/4$ and $2hy=w-u/2+\log3$. Thus $(K_\epsilon)_{ij}=\exp(-\|x_i-x_j\|_2^2)$ exactly, which completes the proof of Theorem 1(b).

The features are abstract points in $\mathbb R^3$, not image patches. As $\epsilon\to0$ the fourth point recedes: its squared distance to the first pair is $-\log(\epsilon/5)$, between $14.5$ and $20$ for $\epsilon$ in the range $10^{-6}$ to $10^{-8}$. The matrix condition numbers stay bounded, but communication between the two groups becomes very weak.

## 8. The mixing gap: a guarantee, an exact classification, and sharpness

### 8.1 A gap-dependent stability bound

**Theorem 2 (gap-dependent step size).** Let $W$ ($n\geq2$) be primitive, stochastic and reversible, with stationary distribution $\pi$ and nonnegative eigenvalues. Let $B$ be symmetric positive semidefinite with $Be=be$, $b>0$. Put $g=1-\lambda_2(W)$ and $M=\sqrt{\kappa_\pi}$. Then both matrices in (1) are Schur stable whenever

\[
0<t\rho(B)<\frac{g}{M},
\tag{27}
\]

and more generally whenever

\[
0<t\rho(B)<\min\left\{2,\frac{g}{M-g}\right\},
\tag{28}
\]

where the fraction is $+\infty$ if $M=g=1$.

Entrywise nonnegativity of $B$ is not needed, and neither is stochasticity. Nonnegative eigenvalues of $W$ hold for every kernel denoiser (2) with $K$ positive semidefinite, including Gaussian kernels, because $W$ is similar to $D^{-1/2}KD^{-1/2}$.

*Proof.* Scaling $B$ by $1/\rho(B)$ and $t$ by $\rho(B)$ leaves (1) unchanged, so assume $\rho(B)=1$. On the quotient modulo constants use the norm

\[
\|[x]\|_\pi=\min_{c\in\mathbb C}\Big(\sum_i\pi_i|x_i-c|^2\Big)^{1/2},\qquad x\in\mathbb C^n.
\]

Let $Q$ and $C$ be the quotients of $W$ and $B$. Reversibility makes $Q$ self-adjoint for this norm, with eigenvalues in $[0,1-g]$. Hence $\|Q\|_\pi=1-g$ and $\|2Q-I\|_\pi\le1$.

*Norm comparison.* For each $c$, $\pi_{\min}\sum_i|x_i-c|^2\le\sum_i\pi_i|x_i-c|^2\le\pi_{\max}\sum_i|x_i-c|^2$. Minimising over $c$ on each side, with possibly different minimisers, gives $\sqrt{\pi_{\min}}\|[x]\|_2\le\|[x]\|_\pi\le\sqrt{\pi_{\max}}\|[x]\|_2$ for the Euclidean quotient norm $\|[x]\|_2$. Hence, for any operator $X$ on the quotient, $\|X\|_\pi\le M\|X\|_2$. The Euclidean quotient of $B$ is symmetric with spectrum in $[0,1]$. Applying the spectral theorem to $\beta$ and to $\beta/(1+t\beta)$, whose maximum on $[0,1]$ is $1/(1+t)$, gives

\[
\|C\|_\pi\leq M,\qquad
\|(I+tC)^{-1}C\|_\pi\leq\frac{M}{1+t}.
\tag{29}
\]

Therefore

\[
\|Q(I-tC)\|_\pi\leq(1-g)(1+tM),
\tag{30}
\]

and, from the identity $R_q-Q=-t(I+tC)^{-1}C(2Q-I)$,

\[
\|R_q\|_\pi\leq1-g+\frac{tM}{1+t}.
\tag{31}
\]

Condition (27) makes both bounds less than one. More generally, (31) is less than one when $t(M-g)<g$, and this also suffices in (30) because $M(1-g)\leq M-g$ for $M\ge1$. The restriction $t<2$ stabilises the constant mode $1-tb$ of $P$; the constant mode $1/(1+tb)$ of $R$ is stable for all $t>0$. $\square$

The technique, a spectral-gap contraction in $L^2(\pi)$ combined with a norm comparison, is standard for reversible Markov chains. We present the bound as a convenient explicit tool. The new point is that its order in $g$ cannot be improved (Section 8.3).

### 8.2 The four-state family: an exact two-parameter classification

We now treat $\epsilon$ and $t$ as independent in (14). This lets us locate exactly where stability is first lost.

The three-dimensional symmetric block of Section 7 contains the constant mode. Removing it with the differences $(x_1-x_4,\,x_3-x_4)$ leaves the PnP-ISTA quotient

\[
N(\epsilon,t)=\frac{Q_0+\epsilon Q_1}{1+\epsilon}(I-tC_4),
\tag{32}
\]

with $Q_0,Q_1,C_4$ given in Appendix A. Its $+1$ condition is governed by

\[
\Delta(\epsilon,t)=\det\big((1+\epsilon)I-(Q_0+\epsilon Q_1)(I-tC_4)\big)
=a_0(\epsilon)+a_1(\epsilon)\,t+a_2(\epsilon)\,t^2,
\tag{33}
\]

an explicit quadratic in $t$ whose coefficients are listed in Appendix B. We use the following elementary test.

**Lemma 1 (quadratic Schur test).** A real monic quadratic $z^2-Tz+D$ has both roots in the open unit disc if and only if $1-D>0$, $1-T+D>0$ and $1+T+D>0$.

*Proof.* If the roots are non-real they are conjugate, $D$ is their squared modulus, and the last two expressions are automatically positive. If they are real, positivity at $z=1$ puts both roots below 1 or both above 1, and the latter would give $D>1$. Positivity at $z=-1$ likewise puts both roots above $-1$. Conversely, suppose both roots lie in the open unit disc. If they are real, they lie in $(-1,1)$, and then $D<1$ and the quadratic is positive at $z=\pm1$, including for repeated roots. If they are non-real, $r$ and $\bar r$, then $D=|r|^2<1$, $1-T+D=|1-r|^2>0$ and $1+T+D=|1+r|^2>0$. This is the degree-two case of the Jury (Schur–Cohn) criterion. $\square$

**Theorem 4 (exact classification of the four-state family).** Let $0\le\epsilon\le10^{-4}$ and $0<t<2$. Then $P=W_\epsilon(I-tB)$ is Schur stable if and only if $\Delta(\epsilon,t)>0$, and the same holds for $R=I-W_\epsilon+(I+tB)^{-1}(2W_\epsilon-I)$. If $\Delta(\epsilon,t)<0$, both have a real eigenvalue greater than one. Moreover:

1. $\Delta(\epsilon,0)>0$ for $\epsilon>0$, and $\Delta(\epsilon,\cdot)$ takes negative values in $(0,2)$ if and only if $\epsilon<\epsilon_{\rm crit}$ (at $\epsilon=\epsilon_{\rm crit}$ it has a double zero, a single neutral step with eigenvalue one), where $\epsilon_{\rm crit}=1.69466328\ldots\times10^{-5}$ is the unique root in $(0,10^{-4}]$ of the discriminant $a_1^2-4a_0a_2$.
2. For $0<\epsilon<\epsilon_{\rm crit}$ the unstable steps in $(0,2)$ form the interval $(T_-(\epsilon),T_+(\epsilon))$ between the two roots of $\Delta(\epsilon,\cdot)$, and $T_+(\epsilon)<2$. The **first loss of stability** is therefore exactly
\[
T_-(\epsilon)=\frac{-a_1-\sqrt{a_1^2-4a_0a_2}}{2a_2}
=\frac{968750}{133}\epsilon+\frac{595515472565625}{4705274}\epsilon^2+O(\epsilon^3),
\tag{34}
\]
and it occurs through the real eigenvalue $+1$, simultaneously for $P$ and $R$.
3. At $\epsilon=0$ (the reducible limit) $\Delta(0,t)=t\,(a_1(0)+a_2(0)t)$ is negative exactly for $0<t<t_0=1330000/2254167\approx0.590$, so both matrices are unstable for every such $t$.
4. Along $\epsilon=t/100000$, the family is unstable exactly for $0<t<t^*$, with $t^*$ as in (22), and Schur stable for $t^*<t<2$.

*Proof.* (i) *Factorisation.* For every $(\epsilon,t)$ the characteristic polynomial of $P$ factors exactly as $(\lambda-(1-t))(\lambda-p_a)(\lambda^2-T_P\lambda+D_P)$. Here $1-t$ is the constant mode, and $p_a=\frac{1/10-3\epsilon/2}{2(1+\epsilon)}(1-t/2500)$ is the antisymmetric mode. The quadratic is the characteristic polynomial of $N(\epsilon,t)$, and $1-T_P+D_P=\Delta/(1+\epsilon)^2$. The analogous factorisation holds for $R$, with constant mode $1/(1+t)$, antisymmetric mode $1-w_a+(2w_a-1)/(1+t/2500)$ where $w_a=\frac{1/10-3\epsilon/2}{2(1+\epsilon)}$, and a quadratic factor $\lambda^2-T_R\lambda+D_R$ with $1-T_R+D_R=\Delta\cdot\theta$ for an explicit rational function $\theta>0$.

(ii) *Side conditions.* On the box $[0,10^{-4}]\times[0,2]$, the four rational functions $1-D_P$, $1+T_P+D_P$, $1-D_R$ and $1+T_R+D_R$, and the factor $\theta$, have numerators and denominators with strictly positive Bernstein coefficients (Appendix C, Table C1). They are therefore positive on the box. The constant modes $1-t$ and $1/(1+t)$ lie in $(-1,1)$ for $0<t<2$, and the antisymmetric modes lie in $(-1,1)$ on the box, again by positive Bernstein coefficients of $1\mp p_a$ and $1\mp r_a$.

(iii) By Lemma 1, $P$ and $R$ are Schur stable exactly when $\Delta>0$. When $\Delta<0$, $\det(I-P)$ and $\det(I-R)$ are negative and the intermediate value argument of Section 7 gives a real eigenvalue greater than one.

(iv) $a_0(\epsilon)=\epsilon(31/120+217\epsilon/200)>0$ and $a_2>0$ on the box. The discriminant has exactly one root in $(0,10^{-4}]$, namely $\epsilon_{\rm crit}$, by the Sturm counts in Table C2. On the box $a_1<0$, so both roots of $\Delta(\epsilon,\cdot)$ are positive, and $4a_2+a_1>0$, so the vertex $-a_1/(2a_2)$ is less than two (Table C1). Moreover $\Delta(\epsilon,2)=a_0+2a_1+4a_2>0$ on the box (Table C1); with the vertex below two, both roots lie in $(0,2)$ when $\epsilon<\epsilon_{\rm crit}$, so $T_+(\epsilon)<2$. This gives items 1 and 2. The expansion in (34) is a direct series computation. Item 3 is $\Delta(0,t)=t\,(a_1(0)+a_2(0)t)$, and item 4 is (19)–(22) combined with the classification.

Steps (i), (ii) and (iv) are exact rational computations performed by `certify_revision.py` and `emit_certificates.py` (Section 15). The decisive polynomials, their Bernstein coefficient arrays, the Sturm counts and the isolating intervals are printed or summarised in Appendix C and stored in hashed coefficient files, which the standard-library checker `check_certificates.py` re-verifies without any symbolic derivation. $\square$

### 8.3 The order in the gap is sharp; the constant is not

The quotient $Q_\epsilon=(Q_0+\epsilon Q_1)/(1+\epsilon)$ has eigenvalues $1-g_\epsilon$ and approximately $47/60$, where

\[
g_\epsilon=\frac{31}{26}\epsilon+O(\epsilon^2),
\qquad
\frac{T_-(\epsilon)}{g_\epsilon}\longrightarrow\frac{812500}{133}\approx6109.
\tag{35}
\]

So no uniform stability guarantee of asymptotically larger order than $g$ is possible, even with a fixed forward matrix, a fixed stationary distribution and the Gaussian conditioning bounds of Theorem 1(b). This sharpness statement concerns $0<\epsilon<\epsilon_{\rm crit}\approx1.7\times10^{-5}$. For larger $\epsilon$ this family has no instability below $t=2$, and the linear approximation (35) is accurate only for $\epsilon\lesssim10^{-6}$. For example, $T_-/\epsilon\approx9171$ at $\epsilon=10^{-5}$.

The constant is far from sharp. Here $\kappa_\pi=2$, so Theorem 2 guarantees $t<g/\sqrt2$, whereas stability actually persists up to about $6109\,g$. The guarantee is conservative by a factor of $\frac{812500}{133}\sqrt2\approx8{,}640$, about 8,600, in the limit. Whether the factor $\sqrt{\kappa_\pi}$ is necessary at all is open. Every unit-diagonal kernel with entries at most one has degrees in $[1,n]$, so $\kappa_\pi\le n$ and $\sqrt{\kappa_\pi}\le\sqrt n$ for the kernels of Theorem 1. The question is therefore only meaningful for large $n$ or for general reversible kernels. We tested it numerically: we optimised over reversible positive kernels in dimensions three and four, with target $\kappa_\pi$ between 5 and 1000. The optimiser found no instance unstable below $t=2$, so the test is inconclusive (Section 14).

### 8.4 The effective negative coefficient

At $\epsilon=0$, the first three coordinates form one component with stationary weights $(5,5,3)/13$, and the fourth is a singleton. Let $U_4$ assign the two component levels to their vertices, and let $\Phi_4$ average the first component with these weights and leave the singleton unchanged. Then $\Phi_4U_4=I$, and the compressed Hessian $\Phi_4BU_4$ has eigenvalues one and

\[
c=\frac1{2500}+\frac{229}{500000}\left(-\frac{16}{13}\right)
=-\frac{133}{812500}.
\tag{36}
\]

The slow mode is damped by about $(31/26)\epsilon$ and amplified by about $|c|\,t$. Equating the two gives exactly the leading term of (34): $(31/26)\epsilon/|c|=\frac{968750}{133}\epsilon$. Along $\epsilon=t/100000$ both slow eigenvalues satisfy

\[
\lambda(t)=1+\frac{1973}{13000000}t+O(t^2).
\tag{37}
\]

## 9. Positive convolution operators

The four-state forward map is not translation invariant. We now impose that restriction, as in deblurring with periodic boundary conditions.

**Theorem 3 (positive circulant forward operator).** There are a fixed strictly positive symmetric circulant stochastic positive definite $A_c\in\mathbb R^{12\times12}$ and genuine Gaussian kernels $K_\epsilon$ such that, with $B_c=A_c^2$ and $\epsilon=t^2$, both operators in (1) have a real eigenvalue greater than one for every

\[
0<t<t_c,\qquad t_c=1.9157074941\ldots\times10^{-4}.
\tag{38}
\]

Here $t_c$ is the unique root in $(0,2]$ of an explicit degree-34 integer polynomial $p(t)$, isolated in $(1.9157074941,\,1.9157074942)\times10^{-4}$ (Appendix C). An exact witness on the curve is $t=10^{-4}$, $\epsilon=10^{-8}$, where $\det(I-P)\approx-2.862\times10^{-17}$ and $\det(I-R)\approx-2.862\times10^{-17}$ (exact fractions in `revision_certificate.json`). For $0<\epsilon\le1/100$ the affinity condition numbers satisfy $\kappa_2(K_\epsilon)<413$ and the stationary-weight ratio is below 4.

The window is narrow and the growth is tiny. For $t_c<t<1$ along the Gaussian curve $\epsilon=t^2$, $\det(I-P)>0$ (the Sturm count of Table C2 excludes any further root of $p$ in $(0,2]$; for $t\ge1$ the algebraic continuation is no longer a Gaussian affinity family), and the numerically computed spectra at $t=5\times10^{-4}$ and $10^{-2}$ are stable. Inside the window the excess spectral radius is about $3.6\times10^{-8}$ (at $t=10^{-4}$). The off-curve finite witness $t=1/100$, $\epsilon=10^{-9}$ (in `circulant_certificate.json`) is also valid, but lies at a much weaker link than the curve $\epsilon=t^2$. On the curve, $t=1/100$ is stable.

### 9.1 A fixed positive circulant matrix

Index coordinates by $0,\ldots,11$, and let $p_i=\cos(\pi i/2)$, $q_i=\sin(\pi i/2)$, $E=ee^{\mathsf T}/12$ and $Z=(pp^{\mathsf T}+qq^{\mathsf T})/6$. Then $p,q$ are orthogonal to each other and to $e$, each with squared norm six, so $Z$ is an orthogonal projection with $EZ=0$ and $Z_{ij}=\cos(\pi(i-j)/2)/6$. Put

\[
A_c=\frac1{12}I+\frac{11}{12}E+\frac5{12}Z,
\tag{39}
\]

a symmetric circulant with first row $\frac1{144}(33,11,1,11,21,11,1,11,21,11,1,11)$. All entries are at least $1/144$ and rows sum to one. Its eigenvalues are $1$, $1/2$ (twice) and $1/12$ (nine times), so $\kappa_2(A_c)=12$ and $\kappa_2(B_c)=144$, where

\[
B_c=A_c^2=\frac1{144}I+\frac{143}{144}E+\frac{35}{144}Z.
\]

No assumption that the convolution weights are Gaussian or decrease with distance is made.

### 9.2 The two affinity components

Partition the indices into $G_1=\{0,1,2,4,6\}$ and $G_2=\{3,5,7,8,9,10,11\}$, with high-similarity subgroups $G_1^h=\{0,2,4,6\}$ and $G_2^h=\{5,8,9,10\}$. Let $K_0$ have unit diagonal, zero entries between components, and within-component entries $99/100$ inside a high-similarity subgroup and $1/100$ otherwise. For $0<\epsilon<1$ set

\[
K_\epsilon=(1-\epsilon)K_0+\epsilon ee^{\mathsf T},\qquad
W_\epsilon=\operatorname{diag}(K_\epsilon e)^{-1}K_\epsilon.
\tag{40}
\]

Each block of $K_0$ is a positive diagonal matrix plus $(98/100)J_h+(1/100)J_{\rm block}$, so $K_0\succeq I/100$ and $K_\epsilon\succeq(1-\epsilon)I/100$. The four degree values at $\epsilon=0$ are $199/50$, $26/25$, $4$ and $53/50$, and the degrees are $(1-\epsilon)d_i+12\epsilon$, so the stationary-weight ratio is at most $50/13<4$. For $0<\epsilon\le1/100$, $\lambda_{\min}(K_\epsilon)\ge0.0099$ and the maximum degree is at most $4.08$, giving $\kappa_2(K_\epsilon)<413$. These bounds do not hold uniformly as $\epsilon\to1$. The stationary vector varies with $\epsilon$ in this construction.

### 9.3 The negative two-group compression

Let $U=[1_{G_1},1_{G_2}]$, and let $\Phi$ average within each group with the stationary weights of $W_0$. Then $\Phi U=I_2$, $W_0U=U$ and $\Phi W_0=\Phi$. On $G_1$ each high vertex has weight $199/848$ and the other $13/212$; on $G_2$ each high vertex has weight $200/959$ and each other $53/959$. Since $Z1_{G_1}=q/6$, whose weighted means on $G_1$ and $G_2$ are $13/1272$ and $241/5754$, and $\Phi EU$ has identical rows, $\Phi B_cU$ has eigenvalues one and

\[
c_{12}=\frac1{144}+\frac{35}{144}\left(\frac{13}{1272}-\frac{241}{5754}\right)
=\frac1{144}-\frac{35}{144}\cdot\frac{12875}{406616}=-\frac{6287}{8364672}<0.
\tag{41}
\]

### 9.4 Proof of Theorem 3

*Asymptotics.* Split $\mathbb R^{12}=\operatorname{ran}U\oplus\ker\Phi$. In this decomposition $W_0$ is block diagonal with blocks $I_2$ and $Q$, where every eigenvalue of $Q$ lies in $(0,1)$. Since $W_{t^2}=W_0+O(t^2)$, a Schur complement gives $\det(I-P(t))=t^2c_{12}\det(I-Q)+O(t^3)$, which is negative for small $t$.

*Certified range.* The asymptotic argument does not say how small $t$ must be, so we compute exactly. Along $\epsilon=t^2$, $\det(I-P(t))=N_c(t)/\det D_{t^2}$, with $N_c(t)=\det(D_{t^2}-K_{t^2}+tK_{t^2}B_c)$, a polynomial of degree 36 with rational coefficients and a factor $t^2$. Let $p(t)$ be $N_c(t)/t^2$ multiplied by a positive integer to clear denominators. Then $p(0)<0$, and the Sturm sign-variation counts of Table C2 show that $p$ has exactly one distinct real root in $(0,2]$, lying in the rational interval $(L,U)$ given there, with $p(L)<0<p(U)$. Define $t_c$ as that root. Hence $p<0$ on $(0,t_c)$. Since $\det D>0$, $\det(I-P)<0$ there. The determinant identity (23), with $B_c$, gives the same sign for $R$. The intermediate value argument then gives real eigenvalues greater than one. $\square$

### 9.5 Exact Gaussian geometry

The off-diagonal affinity levels in (40) are $a_\epsilon=99/100+\epsilon/100$, $b_\epsilon=1/100+99\epsilon/100$ and $\epsilon$, with $1>a_\epsilon>b_\epsilon>\epsilon>0$. Write $\delta_1=-\log a_\epsilon<\delta_2=-\log b_\epsilon<\delta_3=-\log\epsilon$. Take orthonormal vectors $f_0,\ldots,f_{11},u_1,u_2,z$. For a vertex in a high subgroup of $G_k$, set $y_i=\sqrt{\delta_1/2}\,f_i+\sqrt{(\delta_2-\delta_1)/2}\,u_k$; for any other vertex, $y_i=\sqrt{\delta_2/2}\,f_i$. Then

\[
x_i=y_i\pm\frac{\sqrt{\delta_3-\delta_2}}2z\qquad(+\text{ on }G_1,\ -\text{ on }G_2)
\tag{42}
\]

gives squared distances $\delta_1$ within a high subgroup, $\delta_2$ for other pairs in one component, and $\delta_3$ across components. So $(K_\epsilon)_{ij}=\exp(-\|x_i-x_j\|_2^2)$ exactly. The feature dimension here is fifteen; no minimal feature dimension is claimed.

## 10. The unit-diagonal assumption is essential: failure in dimension three

Theorem 1 says that four is the smallest unstable dimension *for unit-diagonal kernels*. This section shows that the restriction is not a technicality. If unequal diagonal entries are allowed, which still gives row-normalised positive kernels, reversible denoisers and a positive forward operator, then instability occurs in dimension three. Dimension three is then minimal (Proposition 1). There is no contradiction with Theorem 1: the identity (7) uses $K_{ii}=1\ge K_{ij}$, and that is exactly what fails here. The family also shows a second phenomenon, stability that is lost and then regained as the step grows.

In this section only, matrices carry a subscript 3.

### 10.1 The operators

Write $E_3=ee^{\mathsf T}/3$ and $v_3=(1,2,-3)^{\mathsf T}$, and fix

\[
B_3=\frac1{180}\begin{pmatrix}68&71&41\\71&86&23\\41&23&116\end{pmatrix}
=\frac{59}{60}E_3+\frac1{60}I+\frac1{30}v_3v_3^{\mathsf T},
\tag{43}
\]

which is strictly positive, symmetric and stochastic, with eigenvalues $1$, $29/60$ and $1/60$. So $\rho(B_3)=1$ and $\kappa_2(B_3)=60$. It comes from a strictly positive forward operator,

\[
A_3=E_3+\frac{I-E_3}{\sqrt{60}}+\frac{\sqrt{29}-1}{14\sqrt{60}}v_3v_3^{\mathsf T},
\tag{44}
\]

with eigenvalues $1$, $\sqrt{29/60}$ and $1/\sqrt{60}$ on the spans of $e$, of $v_3$, and of their orthogonal complement. Hence $A_3^2=B_3$ and $A_3e=e$. Its smallest entry, $(A_3)_{23}=(7\sqrt{60}+2-9\sqrt{29})/(21\sqrt{60})$, is positive, since $\sqrt{60}>15/2$ and $\sqrt{29}<11/2$ make the numerator exceed 5.

Let

\[
\Pi_3=\begin{pmatrix}1&0&0\\0&19/20&1/20\\0&19/20&1/20\end{pmatrix},
\qquad
\pi_3=\frac1{40}(20,19,1)^{\mathsf T},
\]

and, for $0<\epsilon<1$ and $\tau=1-\epsilon$,

\[
W_{3,\epsilon}=\frac{\tau}{2}(I+\Pi_3)+\epsilon e\pi_3^{\mathsf T}.
\tag{45}
\]

$\Pi_3$ and $e\pi_3^{\mathsf T}$ are commuting orthogonal projections in the inner product $\langle x,y\rangle_{\pi_3}=x^{\mathsf T}\operatorname{diag}(\pi_3)y$, with the constants inside the range of $\Pi_3$. So $W_{3,\epsilon}$ has eigenvalues $1$, $\tau$ and $\tau/2$. The term $\epsilon e\pi_3^{\mathsf T}$ makes it strictly positive, and the projection identities give stationarity and detailed balance. Hence

\[
K_{3,\epsilon}=\operatorname{diag}(\pi_3)W_{3,\epsilon}
\tag{46}
\]

is symmetric, positive definite (by congruence with the symmetric matrix $\operatorname{diag}(\pi_3)^{1/2}W_{3,\epsilon}\operatorname{diag}(\pi_3)^{-1/2}$) and entrywise positive, and $W_{3,\epsilon}$ is its row normalisation. Its diagonal is not constant: $(K_{3,\epsilon})_{11}=\frac12-\frac\epsilon4$ while $(K_{3,\epsilon})_{33}=\frac{21-20\epsilon}{1600}$. No Gaussian or unit-diagonal claim is made for this family.

### 10.2 The classification

Stability is controlled by one quadratic,

\[
F_\epsilon(t)
=\frac{\epsilon(1+\epsilon)}2
+\frac{\tau(60\epsilon-1)}{240}t
+\frac{29\tau^2}{7200}t^2.
\tag{47}
\]

**Theorem 5 (three-state reentrant classification).** For $0<\epsilon<1$, with $P_3=W_{3,\epsilon}(I-tB_3)$ and $R_3=I-W_{3,\epsilon}+(I+tB_3)^{-1}(2W_{3,\epsilon}-I)$:

1. for $0<t<2$, $P_3$ is Schur stable if and only if $F_\epsilon(t)>0$;
2. for every $t>0$, $R_3$ is Schur stable if and only if $F_\epsilon(t)>0$;
3. if $F_\epsilon(t)<0$, both have a real eigenvalue greater than one, and if $F_\epsilon(t)=0$, both have eigenvalue one.

The discriminant of $F_\epsilon$ is $\frac{\tau^2}{57600}(3136\epsilon^2-584\epsilon+1)$, which is positive exactly when $\epsilon<\epsilon_3:=\frac{73-\sqrt{3193}}{784}\approx1.76\times10^{-3}$. (The larger root, $\approx0.17$, makes the middle coefficient positive and the roots negative.) Instability in $(0,2)$ therefore occurs only for $\epsilon<\epsilon_3$. For $0<\epsilon\le1/1000$ the two roots

\[
t_\pm(\epsilon)=\frac{15\left(1-60\epsilon\pm\sqrt{3136\epsilon^2-584\epsilon+1}\right)}{29(1-\epsilon)}
\tag{48}
\]

satisfy $0<t_-<t_+<2$. Both iterations are stable on $(0,t_-)\cup(t_+,2)$ and unstable between the roots, and $R_3$ is also stable for all $t\ge2$. Moreover $t_-(\epsilon)=120\epsilon+O(\epsilon^2)$ and $t_+(\epsilon)\to30/29$.

*Proof.* **Contrast coordinates.** In the basis $e=(1,1,1)^{\mathsf T}$, $f_1=(1,-1,-1)^{\mathsf T}$, $f_2=(0,1,-19)^{\mathsf T}$, a vector $x$ has coordinates $\pi_3^{\mathsf T}x$, $x_1-\pi_3^{\mathsf T}x$ and $(x_2-x_3)/20$. In this basis

\[
W_{3,\epsilon}\sim\begin{pmatrix}1&0\\0&\Lambda\end{pmatrix},
\qquad
B_3\sim\begin{pmatrix}1&r^{\mathsf T}\\0&C_3\end{pmatrix},
\qquad
\Lambda=\begin{pmatrix}\tau&0\\0&\tau/2\end{pmatrix},\quad
C_3=\begin{pmatrix}-1/120&-59/80\\1/60&61/120\end{pmatrix},
\tag{49}
\]

with $r^{\mathsf T}=(-17/72,\,-767/240)$, which does not affect eigenvalues. Although $C_3$ has positive eigenvalues, $(C_3)_{11}=-1/120<0$; this is the oblique compression of Section 4 in its simplest form. Write $\sigma_1=\operatorname{tr}\Lambda=\tfrac32\tau$, $\sigma_2=\det\Lambda=\tfrac12\tau^2$, $\gamma_1=\operatorname{tr}C_3=\tfrac12$, $\gamma_2=\det C_3=\tfrac{29}{3600}$ and $\mu=\operatorname{tr}(\Lambda C_3)=\tfrac{59}{240}\tau$.

**$P_3$.** The contrast matrix $\Lambda(I-tC_3)$ has trace $T_P=\sigma_1-t\mu$ and determinant $D_P=\sigma_2(1-\gamma_1t+\gamma_2t^2)$, and $1-T_P+D_P=F_\epsilon(t)$. For $0<t<2$, $T_P=\tau(\tfrac32-\tfrac{59}{240}t)>0$ and $0<D_P=\sigma_2(1-\tfrac t{60})(1-\tfrac{29t}{60})<1$, so the other two conditions of Lemma 1 hold. The constant eigenvalue $1-t$ is stable for $0<t<2$. This proves item 1, and also $\det(I-P_3)=tF_\epsilon(t)$.

**$R_3$.** With $G(t)=\det(I+tC_3)=1+\gamma_1t+\gamma_2t^2>0$, the contrast matrix $I-\Lambda+(I+tC_3)^{-1}(2\Lambda-I)$ has trace $T_N/G$ and determinant $D_N/G$, where

\[
T_N=\sigma_1+t\{\gamma_1(1+\sigma_1)-2\mu\}+\gamma_2t^2(2-\sigma_1),\qquad
D_N=\sigma_2+t\{\gamma_1(\sigma_1-\sigma_2)-\mu\}+\gamma_2t^2(1-\sigma_1+\sigma_2),
\]

using $(I+tC_3)^{-1}=((1+t\gamma_1)I-tC_3)/G$ and $C_3^2=\gamma_1C_3-\gamma_2I$. Then $G-T_N+D_N=F_\epsilon(t)$. The other two conditions,

\[
G-D_N=1-\sigma_2+t\{\gamma_1(1-\sigma_1+\sigma_2)+\mu\}+\gamma_2t^2(\sigma_1-\sigma_2)
\]

and

\[
G+T_N+D_N=1+\sigma_1+\sigma_2+t\{\gamma_1(2+2\sigma_1-\sigma_2)-3\mu\}+\gamma_2t^2(4-2\sigma_1+\sigma_2),
\]

are positive for all $t>0$, because $\sigma_2<1$, $1-\sigma_1+\sigma_2=(1-\tau)(1-\tau/2)>0$, $\sigma_1-\sigma_2=\tau(3-\tau)/2>0$, $\gamma_1(2+2\sigma_1-\sigma_2)-3\mu=1+\tfrac{61}{80}\tau-\tfrac14\tau^2>0$ and $4-2\sigma_1+\sigma_2=4-3\tau+\tfrac12\tau^2>0$. The constant eigenvalue $1/(1+t)$ is stable. This proves item 2.

**Item 3** follows from $\det(I-R_3)=\det(I-P_3)/\det(I+tB_3)=tF_\epsilon(t)/((1+t)G(t))$ and the intermediate value argument of Section 7.

**Roots.** For $0<\epsilon\le1/1000$, $3136\epsilon^2-584\epsilon+1\ge1-0.584>0$, $(1-60\epsilon)^2-(3136\epsilon^2-584\epsilon+1)=464\epsilon(1+\epsilon)>0$ and $1-60\epsilon>0$, so both roots are positive. Moreover $t_+<30/(29(1-\epsilon))\le30000/28971<2$. Since the leading coefficient of $F_\epsilon$ is positive, $F_\epsilon<0$ exactly between the roots. Finally, $\partial_\epsilon F(0,0)=1/2$ and $\partial_tF(0,0)=-1/240$, so $t_-'(0)=120$. $\square$

**A concrete witness.** For $\epsilon=1/1000$ and $t=1/4$, $\det(I-P_3)=-2898619/51200000000$ and $\det(I-R_3)=-26087571/648290000000$. The thresholds are $t_-\approx0.15149$ and $t_+\approx0.82189$.

**What the quantifiers imply.** The forward operator, $B_3$ and $\pi_3$ are fixed, and for $\epsilon\le1/1000$ every denoiser eigenvalue is at least $999/2000$. Yet $g_\epsilon=\epsilon$ and $t_-(\epsilon)=120g_\epsilon+O(g_\epsilon^2)$. So no positive step guarantee can depend only on the dimension, $B$, $\pi$ and a lower bound on the denoiser's eigenvalues; the mixing gap must enter.

**Proposition 1 (dimension three is minimal without unit diagonals).** Let $n\le2$, let $W$ be primitive and stochastic, and let $B$ be symmetric positive semidefinite and stochastic. Then both operators are Schur stable for $0<t<2$.

*Proof.* Dimension one is scalar. For $n=2$, $B=\begin{pmatrix}a&1-a\\1-a&a\end{pmatrix}$ has eigenvalues $1$ and $\beta=2a-1\in[0,1]$, and on the quotient $W$ acts by its other eigenvalue $w\in(-1,1)$. The nonconstant eigenvalues of $P$ and $R$ are $w(1-z)$ and $(w+z(1-w))/(1+z)$ with $z=t\beta\in[0,2)$. The first has modulus less than one. For the second, $1-\frac{w+z(1-w)}{1+z}=\frac{1-w(1-z)}{1+z}>0$ and $1+\frac{w+z(1-w)}{1+z}=\frac{1+w+z(2-w)}{1+z}>0$. $\square$

The case of Proposition 1 with positive $W$ is also covered by Banerjee and Chaudhury (2025, Proposition 1.2), since every symmetric stochastic $2\times2$ matrix has the form $\alpha I+\beta ee^{\mathsf T}$. We include the proof only to make the minimality statement self-contained.

## 11. Symmetric denoisers: stability in every dimension, with one obstruction

The obstruction in Sections 4–10 depends on the denoiser being non-symmetric. When the denoiser is symmetric, the oblique compression becomes orthogonal, and in any dimension the only possible failure is a vector fixed by $W$ and annihilated by $B$. This holds even if the denoiser is not stochastic, as for the symmetric normalisation $D^{-1/2}KD^{-1/2}$.

**Theorem 6 (symmetric denoisers).** Let $W\in\mathbb R^{n\times n}$ be symmetric with spectrum in $(-1,1]$. Let $B$ be symmetric positive semidefinite, and let $0<t<2/\rho(B)$. Then $P(t)$ and $R(t)$ in (1) are Schur stable if and only if there is no nonzero $x$ with $Wx=x$ and $Bx=0$.

*Proof.* If $Wx=x$ and $Bx=0$ with $x\ne0$, then $P(t)x=x$ and $R(t)x=x-x+(2x-x)=x$, so neither is Schur stable.

Conversely, suppose no such $x$ exists. All norms are Euclidean on $\mathbb C^n$, and since $W$ is symmetric, $\|W\|_2=\rho(W)\le1$. Eigenvectors of $B$ with eigenvalue $\beta$ satisfy $|1-t\beta|<1$ unless $\beta=0$, because $0<t\beta<2$ when $\beta>0$.

*$P$.* Let $P(t)x=\lambda x$ with $x\ne0$ and $|\lambda|\ge1$. Then $\|x\|\le|\lambda|\|x\|=\|W(I-tB)x\|\le\|(I-tB)x\|\le\|x\|$, so equality holds throughout. Equality in the last step forces every component of $x$ along eigenvectors of $B$ with $\beta>0$ to vanish, so $Bx=0$ and $(I-tB)x=x$. Then $Wx=\lambda x$ with $|\lambda|\ge1$ and $\lambda\in\sigma(W)\subset(-1,1]$, so $\lambda=1$, a contradiction.

*$R$.* Let $R(t)x=\lambda x$ with $x\ne0$ and $|\lambda|\ge1$. Multiplying by $I+tB$ gives (9) with $B$ in place of $C$ and $W$ in place of $Q$: $[\lambda I+(\lambda-1)tB]x=(I-tB)Wx$. As in Section 5.2, the left-hand matrix is invertible and $x=\varphi_\lambda(B)Wx$ with $\|\varphi_\lambda(B)\|_2\le1$. Hence $\|x\|\le\|Wx\|\le\|x\|$. Since $W$ is symmetric, $\|Wx\|=\|x\|$ forces $x$ to lie in the span of eigenvectors with eigenvalues of modulus one, and since $-1\notin\sigma(W)$, $Wx=x$. Then $x=\varphi_\lambda(B)x$, so $\varphi_\lambda(\beta)=1$ for every eigenvalue $\beta$ of $B$ carrying a nonzero component of $x$. That means $\lambda(1+t\beta)=1$, which with $|\lambda|\ge1$ forces $\beta=0$ and $\lambda=1$. Thus $Bx=0$ and $Wx=x$, a contradiction. $\square$

**Corollary 1.** Let $n\ge2$ and let $K$ be symmetric, nonnegative and irreducible (so $D=\operatorname{diag}(Ke)$ is positive), and suppose either $K\succeq0$ or $D^{-1}K$ is primitive. Let $B\succeq0$ with $Be=be$ and $b>0$. Then the symmetrically normalised denoiser $W_s=D^{-1/2}KD^{-1/2}$ gives Schur stable $P(t)$ and $R(t)$ for all $0<t<2/\rho(B)$, in every dimension.

*Proof.* $W_s$ is similar to $D^{-1}K$, so its spectrum lies in $[-1,1]$. The eigenvalue $-1$ is excluded by positive semidefiniteness or by primitivity. The eigenvalue $1$ is simple, with eigenvector $D^{1/2}e$, which has positive entries. If $Bx=0$ then $b\,e^{\mathsf T}x=(Be)^{\mathsf T}x=e^{\mathsf T}Bx=0$, so $\ker B\perp e$. But $e^{\mathsf T}D^{1/2}e>0$, so $D^{1/2}e\notin\ker B$, and Theorem 6 applies. $\square$

**Credit.** Theorem 6 is close to known results. Athalye, Chaudhury and Kumar (2023) proved, for symmetric stochastic positive semidefinite $W$, that $\|P\|_2<1$ and $\|R\|_2<1$ if and only if $\ker(I-W)\cap\ker A=\{0\}$, where $\ker A=\ker B$ (their Theorems 1 and 2), via a lemma on products of symmetric nonexpansive matrices that does not itself need stochasticity. Banerjee and Chaudhury's Theorem 1.2 is adapted from that work, and our spectral-norm argument follows it. Athalye et al. also used the symmetric normalisation $D^{-1/2}KD^{-1/2}$ as a proof device for inpainting. What Theorem 6 adds is modest: the same characterisation, stated for Schur stability, when $W$ is symmetric but need not be stochastic or positive semidefinite (spectrum in $(-1,1]$), including the PnP-ADMM operator; and Corollary 1, which records that the symmetrically normalised denoiser of any irreducible kernel is stable on the full step interval when $Be=be$. As a numerical check, the symmetric normalisation of the four-state kernel at $t=1/4$ has spectral radius $1-1.03\times10^{-4}$, whereas the row-normalised denoiser has $1+1.71\times10^{-5}$. Note that $W_s$ does not preserve constants, so it changes the denoiser and its fixed point; we recommend it only as a stability safeguard, alongside the doubly stochastic symmetrisation of Sreehari et al. (2016) and the scaled updates of Sinha et al. (2025).

## 12. What the results mean in practice

**The instability is real but extremely slow.** The excess spectral radius of the four-state family is at most $1.72\times10^{-5}$ per iteration (Figure 1(b)). In 50-digit arithmetic:

| step $t$ | $\rho(P)-1$ | iterations for tenfold growth | smallest affinity |
|---|---|---|---|
| $10^{-2}$ | $1.48\times10^{-6}$ | $1.6\times10^6$ | $2\times10^{-8}$ |
| $0.1$ | $1.15\times10^{-5}$ | $2.0\times10^5$ | $2\times10^{-7}$ |
| $0.25$ | $1.71\times10^{-5}$ | $1.3\times10^5$ | $5\times10^{-7}$ |
| $0.5$ | $3.32\times10^{-6}$ | $6.9\times10^5$ | $1\times10^{-6}$ |

$R$ agrees with $P$ to the digits shown. Tenfold growth needs roughly $10^5$ to $10^6$ iterations, against the 20 to 500 iterations of a typical PnP run. In the twelve-state convolution example the excess is about $3.6\times10^{-8}$. Within the iteration counts actually used, the unstable direction would be invisible, and it is below single-precision rounding. The results are structural: they show that no *uniform* guarantee of the conjectured form exists. They do not demonstrate divergence in an imaging pipeline.

**The affinities involved are very small.** The weak link has affinity between $10^{-8}$ and $10^{-6}$ in the four-state family, and below about $4\times10^{-8}$ in the convolution family. Practical non-local-means and bilateral implementations use $k$-nearest-neighbour lists, search windows or thresholds that set affinities this small to zero. Truncation does not remove the danger, it moves it. A truncated kernel can be reducible, and then $g=0$, Banerjee and Chaudhury's Theorem 1.1 no longer applies, and Theorem 4 at $\epsilon=0$ shows that the reducible limit is unstable for all $0<t<0.59$.

**The features are abstract points.** The kernels are built from feature vectors in $\mathbb R^3$ (four-state) and $\mathbb R^{15}$ (twelve-state), not from image patches, and no spatial search window is imposed. Whether patch-based kernels on natural images can realise the same geometry is open.

**Gap-aware step-size advice.** Practitioners using row-normalised kernel denoisers can proceed as follows.

1. Compute or estimate the mixing gap $g=1-\lambda_2(W)$ and the weight imbalance $\kappa_\pi=\max_id_i/\min_id_i$, where $d=Ke$. The gap is the second eigenvalue of the symmetric matrix $D^{-1/2}KD^{-1/2}$, available from a few Lanczos iterations.
2. If $g$ is not small, the iteration is safe for $t\rho(B)<\min\{2,g/(\sqrt{\kappa_\pi}-g)\}$ (Theorem 2). This bound is conservative: in our example by a factor of about 8,600.
3. If $g$ is very small, which happens with nearly disconnected similarity graphs, sharply separated textures, small bandwidths or aggressive truncation, either reduce the step in proportion to $g$, or remove the obstruction. Symmetrising the denoiser makes every $t<2/\rho(B)$ safe in every dimension (Theorem 6, Corollary 1). Scaled PnP updates (Sinha et al., 2025) achieve the same with linear convergence.
4. A very small $g$ also means slow convergence of the averaging itself, whose contraction factor is about $1-g$, so a small gap is worth detecting regardless of stability.

**Mapping the assumptions onto imaging problems.**

- *Deblurring with periodic boundaries.* A blur with nonnegative weights summing to one is a circulant stochastic $A$, and $B=A^{\mathsf T}A$ satisfies $Be=e$ and $\rho(B)=1$. Theorems 1(a), 2, 3 and 6 all apply. Theorem 3 shows that a positive blur alone does not prevent instability.
- *Cartesian MRI.* With unitary DFT $F$ and sampling mask $S$, $B=F^{\mathsf H}SF$ is a circulant orthogonal projection. Since $Fe$ is supported on the zero frequency, $Be=e$ holds exactly when the DC sample is acquired, which is standard. Then $b=\rho(B)=1$. $B$ is real symmetric when the mask is conjugate-symmetric, and otherwise is complex Hermitian, which is outside our real setting. Theorems 1(a), 2 and 6 apply. $B$ is not entrywise positive, so our emphasis on positive forward operators matters less here. The Hessian class of Theorem 1(a) includes such projections.
- *Inpainting.* $B$ is a diagonal 0/1 mask, so $Be\le e$ but $Be\ne be$ in general. Theorem 1(a) and Theorem 2 do not apply. However, Banerjee and Chaudhury's Proposition 1.1 already shows stability on the full interval for every irreducible stochastic $W$ when $B$ is diagonal, in every dimension. Our counterexamples, which have non-diagonal $B$, do not touch inpainting.
- *Learned denoisers.* CNN denoisers such as DRUNet are nonlinear and not row-stochastic. Our results bear on them only through linearisation at a fixed point, and we make no claim about them.

## 13. Scope, limitations and non-claims

All qualifications are collected here.

- **Quantifiers.** In Theorems 1(b) and 3 the kernel is chosen for each step and then held fixed. Every fixed primitive kernel is stable for all sufficiently small steps (Banerjee and Chaudhury, 2025, Theorem 1.1). Theorems 4 and 5 describe the exact behaviour of each fixed kernel in their families.
- **Class.** "Sharp dimension four" refers to irreducible unit-diagonal kernels with entries at most one and Hessians with $Be=be$. It fails for reducible kernels (Section 3) and for non-unit diagonals (Section 10). It is untested for $Be\le\rho(B)e$.
- **Range.** The four-dimensional instability is certified for $0<t<t^*\approx0.5313$ in our family. Larger steps in dimension four are open. The convolution instability is certified only for $0<t<1.92\times10^{-4}$ along $\epsilon=t^2$.
- **Sharpness of Theorem 2.** The order $g$ is sharp, but only within one family and for $\epsilon<1.7\times10^{-5}$. The constant is off by about 8,600, and the role of $\sqrt{\kappa_\pi}$ is untested.
- **Practical relevance.** The growth rates are $10^{-5}$ to $10^{-8}$ per step, the affinities are tiny, and the features are abstract points. No imaging experiment is reported. See Section 12.
- **Not the first disproof.** Patel's post gave a counterexample to Conjecture 1.1 first; MathDB labels it unverified, and a referee of this paper checked its determinant identities. Gavaskar et al. (2021) exhibited ADMM divergence with a signed forward model.
- **Symmetric denoisers.** Theorem 6 extends Athalye, Chaudhury and Kumar (2023); its stability conclusion needs the condition $\ker(I-W)\cap\ker B=\{0\}$, and without it symmetric denoisers can fail (for example $W=I_2$, $B=ee^{\mathsf T}/2$).
- **Prior art.** Our literature search was bounded, through 27 October 2026. Liu, Xing and Sun (2021) was read in abstract only. Citation-forward searches on Banerjee and Chaudhury (2025) returned no competing resolution, but that is weak evidence. Sinha et al.'s open-question sentence was checked in the arXiv version (2505.15318v1), not in the published version. Novelty claims are to the best of our knowledge.
- **Verification.** The proofs are written out in full. The large-integer identities, the Bernstein sign certificates and the Sturm root counts are exact rational computations by the scripts in Section 15; their decisive data are in Appendix C and in hashed coefficient files. They are not formally verified in a proof assistant. Floating-point numbers in Section 12 and Figure 1 are illustrations, not proofs.

## 14. Open problems

1. **Large steps in dimension four.** Is there an irreducible unit-diagonal (Gaussian) kernel on four points, with $B\succeq0$ and $Be=e$, that is unstable at some step $t\in(t^*,2)$? Our exploratory numerical search found apparent instabilities in this range, but only for degenerate kernels: some affinities were zero or below $10^{-150}$, and $B$ was singular. These are not counterexamples in the class of Theorem 1, and no conclusion is drawn.
2. **Is $\sqrt{\kappa_\pi}$ needed in Theorem 2?** For unit-diagonal kernels $\kappa_\pi\le n$. For general reversible kernels the question is whether the first unstable step can be $o(g)$ as $\kappa_\pi\to\infty$. Our numerical optimisation in dimensions three and four found no unstable instance, which is inconclusive.
3. **The weaker hypothesis $Be\le\rho(B)e$** of Conjecture 1.1, as in super-resolution: does Theorem 1(a) survive in dimension three?
4. **The minimal dimension for positive circulant forward operators** (currently at most twelve), and a circulant family with a wider certified window.
5. **A single fixed Gaussian kernel unstable on a step interval.** Theorem 4 provides one: any fixed $0<\epsilon<\epsilon_{\rm crit}$. Can the interval be made long, and the growth large?
6. **Patch realisations and worst-case growth.** Can the obstruction be realised with natural-image patches and a search window? How large can $\rho-1$ be at a given mixing gap?

## 15. Reproducibility

Each script uses exact rational arithmetic unless stated otherwise. Run them from a copy of the package; each writes its outputs beside itself (or to `--out`).

~~~text
uv run python3 -I verify_exact.py         # Thm 1(b)
uv run python3 -I verify_circulant.py     # Thm 3 constants
uv run python3 -I verify_reentrant.py     # Thm 5
uv run --with sympy==1.13.3 \
    python3 -I certify_revision.py        # Thms 3, 4; (16), (22)
uv run --with sympy==1.13.3 \
    python3 -I practical_numbers.py       # Section 12 table
uv run --with numpy==2.1.3 \
    --with matplotlib==3.9.2 \
    python3 -I make_figure.py             # Figure 1
uv run --with sympy==1.13.3 \
    python3 -I emit_certificates.py       # Appendix C files
uv run python3 -I check_certificates.py   # Appendix C check
~~~

The first three, and `check_certificates.py`, need only the Python standard library. Section headers in `verification_output.txt` separate the output of each script. Tested with Python 3.11.

## Declarations

**Authorship.** Anonymous.

**AI provenance.** This manuscript, its proofs, its verification scripts and its figure were produced with substantial assistance from AI language-model agents. The agents drafted the text, constructed the examples, wrote and ran the exact-arithmetic verifiers, and revised the paper in response to five simulated referee reports and an independent automated proof check. A human editor set the scope and reviewed the outputs. No independent human peer review has yet taken place. Readers should rely on the exact certificates and written proofs rather than on the provenance.

**Competing interests.** None.

**Data and code availability.** All scripts and their outputs accompany this manuscript (Section 15). They are archived at https://doi.org/10.5281/zenodo.23690125 and available at https://github.com/ipitchford/gaussian-pnp-stability (tag `v1.0.0-candidate`).

## Appendix A. Small matrices for checking the algebra

**A.1 The three-dimensional invariant block.** On vectors with $x_1=x_2$, combine the first two columns and use rows $1,3,4$:

\[
\bar K_0=\begin{pmatrix}19/10&1/10&0\\1/5&1&0\\0&0&1\end{pmatrix},
\qquad
\bar H=\begin{pmatrix}3/2&3/10&1/5\\3/5&0&3/5\\2/5&3/5&0\end{pmatrix},
\]

\[
\bar B=\begin{pmatrix}
69853/125000&46903/250000&63391/250000\\
46903/125000&191331/500000&121057/500000\\
63391/125000&121057/500000&125379/500000
\end{pmatrix},
\]

$\bar D_t=(1+\epsilon)\operatorname{diag}(2,6/5,1)$ and $\bar M_t=\bar D_t-(\bar K_0+\epsilon\bar H)+t(\bar K_0+\epsilon\bar H)\bar B$. Substituting $\epsilon=t/100000$ gives (19)–(20).

**A.2 The two-dimensional quotient.**

\[
Q_0=\begin{pmatrix}19/20&1/20\\1/6&5/6\end{pmatrix},\qquad
Q_1=\begin{pmatrix}7/20&-9/20\\1/10&-3/5\end{pmatrix},\qquad
C_4=\begin{pmatrix}
6462/125000&-27251/500000\\
-16488/125000&35137/250000
\end{pmatrix}.
\]

## Appendix B. Exact constants

The cubic in (20) is

\[
f(t)=-\frac{1973}{25000000}
+\frac{74272809}{500000000000}t
-\frac{236896251}{625000000000000000}t^2
-\frac{4748139}{1562500000000000000000}t^3.
\]

The coefficients of $\Delta$ in (33) are

\[
\begin{aligned}
a_0&=\epsilon\left(\frac{31}{120}+\frac{217}{200}\epsilon\right),\qquad
a_1=-\frac{133}{3750000}+\frac{1782779}{10000000}\epsilon+\frac{969583}{50000000}\epsilon^2,\\
a_2&=\frac{751389}{12500000000}-\frac{15987}{1000000000}\epsilon-\frac{1582713}{125000000000}\epsilon^2.
\end{aligned}
\]

The linear terms of $\Delta$, together with $1-47/60=13/60$, give $g_\epsilon$ and the leading term of $T_-$ in (34)–(35).

## Appendix C. Certificates for Theorems 3 and 4

This appendix lets a reader check Theorems 3 and 4 from the paper and three small data files, without rerunning the symbolic derivation. `emit_certificates.py` (sympy) writes the files. `check_certificates.py` uses only the Python standard library. It reads the files, recomputes every Bernstein array and Sturm chain in exact rational arithmetic, and prints the hashes below.

- `certificates/thm4_classification.json`: $T_P,D_P,T_R,D_R,\theta$; numerators, denominators and Bernstein arrays of every side condition; $\Delta$; the discriminant. SHA-256 `f783ed0bc72e689a`\allowbreak`3734599defb009cd`\allowbreak`662b056d392a5aa4`\allowbreak`33cdff9f759ff5ba`.
- `certificates/thm4_tstar.json`: the integer cubic for $t^*$ and its isolating interval. SHA-256 `f00414decf125afb`\allowbreak`08a335eaae23b69c`\allowbreak`f6dc36bf40fac208`\allowbreak`7458dfd2d3f181f6`.
- `certificates/thm3_circulant.json`: the integer coefficients of $p(t)=N_c(t)/t^2$ (degree 34) and the isolating interval for $t_c$. SHA-256 `4badd51b36eee7a9`\allowbreak`34f5f02b50f6ede6`\allowbreak`3b0f62731c8a23c5`\allowbreak`b31e27b6f8428d30`.

**C.1 The quadratic factors of Theorem 4.** With $\epsilon$ and $t$ as in Section 8.2,

\[
\begin{aligned}
T_P&=\frac{184929\,\epsilon t-3750000\,\epsilon-2258335\,t+26750000}{15000000\,(1+\epsilon)},\\
D_P&=\frac{(2500-t)(47961t-250000)(99\epsilon^2+125\epsilon-470)}{375000000000\,(1+\epsilon)^2},\\
\theta&=\frac{625000000}{(1+\epsilon)^2\,(t+2500)\,(47961t+250000)},
\end{aligned}
\]

and $1-T_R+D_R=\theta\,\Delta$. The longer expressions for $T_R$ and $D_R$ have denominators $60(1+\epsilon)(t+2500)(47961t+250000)$ and $600(1+\epsilon)^2(t+2500)(47961t+250000)$. They are given exactly in the first file.

**C.2 Bernstein sign summaries (Table C1).** Each condition is written as a numerator over a denominator, with signs chosen so that the denominator is positive. Substituting $\epsilon=u/10000$ and $t=2s$ maps the box $[0,10^{-4}]\times[0,2]$ onto $[0,1]^2$. Each polynomial is then expanded in the tensor Bernstein basis of its own bidegree. A single box suffices for every condition, and no subdivision is needed. Strictly positive coefficients imply positivity on the closed box.

| Condition | Numerator: coefficients, minimum | Denominator: coefficients, minimum |
|---|---|---|
| $1-D_P$ | 9, $8.13\times10^{10}$ | 3, $3.75\times10^{11}$ |
| $1+T_P+D_P$ | 9, $1.11\times10^{12}$ | 3, $3.75\times10^{11}$ |
| $1-D_R$ | 9, $8.13\times10^{10}$ | 9, $3.75\times10^{11}$ |
| $1+T_R+D_R$ | 9, $1.34\times10^{12}$ | 9, $3.75\times10^{11}$ |
| $\theta$ | 1, $6.25\times10^{8}$ | 9, $6.25\times10^{8}$ |
| $1-p_a$, $1+p_a$ | 4, $4.75\times10^{4}$; 4, $5.25\times10^{4}$ | 2, $5.00\times10^{4}$ |
| $1-r_a$, $1+r_a$ | 4, $4.75\times10^{4}$; 4, $5.25\times10^{4}$ | 4, $5.00\times10^{4}$ |
| $a_2$ | 3, $7.51\times10^{6}$ | 1, $1.25\times10^{11}$ |
| $-a_1$ | 3, $2.65\times10^{3}$ | 1, $1.50\times10^{8}$ |
| $4a_2+a_1$ (vertex below 2) | 3, $9.61\times10^{6}$ | 1, $4.69\times10^{10}$ |
| $\Delta(\epsilon,2)=a_0+2a_1+4a_2$ | 3, $1.59\times10^{7}$ | 1, $9.38\times10^{10}$ |

The minima are of the integer-scaled polynomials stored in the file, so only their signs matter.

**C.3 Sturm counts and isolating intervals (Table C2).** The three certified roots are listed below. $V(x)$ is the number of sign variations of the Sturm chain at $x$. The number of distinct real roots in $(a,b]$ is $V(a)-V(b)$.

- *Discriminant* $a_1^2-4a_0a_2$, cleared to $28302400-1682115068400\,\epsilon+709588614431865\,\epsilon^2+157425214172130\,\epsilon^3+9697236140601\,\epsilon^4$. Isolating interval $\epsilon_{\rm crit}\in(L,U)=(1394/82258229,\ 343/20240009)$, both endpoints $\approx1.694663\times10^{-5}$. Counts $V(0)=3$, $V(L)=3$, $V(U)=2$, $V(10^{-4})=2$: exactly one root in $(0,10^{-4}]$; the polynomial is square-free and changes sign across $(L,U)$.
- *The $t^*$ cubic* $\Delta(t/10^5,t)/t$, cleared to $-123312500000000000+232102528125000000\,t-592240627500\,t^2-4748139\,t^3$. Isolating interval $t^*\in(0.5312853253,\ 0.5312853254)$. Counts $V(0)=2$, $V(L)=2$, $V(U)=1$, $V(2)=1$: negative on $(0,t^*)$, positive on $(t^*,2]$.
- *The circulant polynomial* $p(t)=N_c(t)/t^2$, cleared, degree 34, $p(0)<0$. Isolating interval $t_c\in(1.9157074941,\ 1.9157074942)\times10^{-4}$. Counts $V(0)=19$, $V(L)=19$, $V(U)=18$, $V(1)=18$, $V(2)=18$: the only distinct root in $(0,2]$, with $p(L)<0<p(U)$.

The degree-34 polynomial has repeated factors: $\gcd(p,p')$ has degree 6. The chain therefore ends at that gcd and counts distinct roots. This does not affect the conclusion, because no root at all lies in $(0,L]$, and $p$ changes sign across $(L,U)$. The exact witness at $t=10^{-4}$, $\epsilon=10^{-8}$ is given as two exact fractions in `revision_certificate.json`.

**C.4 What the checker does not do.** It takes the stored polynomials as given. That they are the polynomials of the matrices in Sections 6–9 is established by `certify_revision.py` and `emit_certificates.py`, which build them from the matrix definitions. An independent rebuild from the matrix definitions, written for the referee's proof audit, reproduced the isolating intervals for $t^*$ and $t_c$ and the value $p(0)<0$.

## References

Athalye, C. D., Chaudhury, K. N., & Kumar, B. (2023). On the contractivity of plug-and-play operators. *IEEE Signal Processing Letters, 30*, 1447–1451. https://doi.org/10.1109/LSP.2023.3323248

Banerjee, D., & Chaudhury, K. N. (2025). *On the Schur stability of some image reconstruction operators*. arXiv:2509.06715. https://arxiv.org/abs/2509.06715

Chan, S. H., Wang, X., & Elgendy, O. A. (2017). Plug-and-play ADMM for image restoration: Fixed-point convergence and applications. *IEEE Transactions on Computational Imaging, 3*(1), 84–98. https://doi.org/10.1109/TCI.2016.2629286

Gavaskar, R. G., & Chaudhury, K. N. (2020). Plug-and-play ISTA converges with kernel denoisers. *IEEE Signal Processing Letters, 27*, 610–614. https://doi.org/10.1109/LSP.2020.2986643

Gavaskar, R. G., Athalye, C. D., & Chaudhury, K. N. (2021). On plug-and-play regularization using linear denoisers. *IEEE Transactions on Image Processing, 30*, 4802–4813. https://doi.org/10.1109/TIP.2021.3075092

Hurault, S., Leclaire, A., & Papadakis, N. (2021). *Gradient step denoiser for convergent plug-and-play*. arXiv:2110.03220. https://arxiv.org/abs/2110.03220

Hurault, S., Leclaire, A., & Papadakis, N. (2022). *Proximal denoiser for convergent plug-and-play optimization with nonconvex regularization*. arXiv:2201.13256. https://arxiv.org/abs/2201.13256

Liu, T., Xing, L., & Sun, Z. (2021). Study on convergence of plug-and-play ISTA with adaptive-kernel denoisers. *IEEE Signal Processing Letters, 28*, 1918–1922. https://doi.org/10.1109/LSP.2021.3111594 [abstract only consulted]

Pesquet, J.-C., Repetti, A., Terris, M., & Wiaux, Y. (2021). Learning maximally monotone operators for image recovery. *SIAM Journal on Imaging Sciences, 14*(3), 1206–1237. https://doi.org/10.1137/20M1387961

Romano, Y., Elad, M., & Milanfar, P. (2017). The little engine that could: Regularization by denoising (RED). *SIAM Journal on Imaging Sciences, 10*(4), 1804–1844. https://doi.org/10.1137/16M1102884

Ryu, E. K., Liu, J., Wang, S., Chen, X., Wang, Z., & Yin, W. (2019). *Plug-and-play methods provably converge with properly trained denoisers*. arXiv:1905.05406. https://arxiv.org/abs/1905.05406

Sinha, A., Kumar, B., Athalye, C. D., & Chaudhury, K. N. (2025). Linear convergence of plug-and-play algorithms with kernel denoisers. *IEEE Transactions on Signal Processing, 73*, 2646–2659. https://doi.org/10.1109/TSP.2025.3573044 (arXiv:2505.15318)

Sreehari, S., Venkatakrishnan, S. V., Wohlberg, B., Buzzard, G. T., Drummy, L. F., Simmons, J. P., & Bouman, C. A. (2016). Plug-and-play priors for bright field electron tomography and sparse interpolation. *IEEE Transactions on Computational Imaging, 2*(4), 408–423. https://doi.org/10.1109/TCI.2016.2599778
