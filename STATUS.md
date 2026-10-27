# Status

Version 1.0.0-candidate (27 October 2026). DOI [10.5281/zenodo.23690125](https://doi.org/10.5281/zenodo.23690125) (concept DOI [10.5281/zenodo.23690124](https://doi.org/10.5281/zenodo.23690124)). Scholarly creator: Anonymous.

**This is an unrefereed research candidate.** No journal and no named human specialist has reviewed it. The only checking was internal and AI-assisted (see [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md)). That is not peer review, independent reproduction or formal verification.

**Target.** The open question in Sinha, Kumar, Athalye and Chaudhury, *IEEE Trans. Signal Process.* 73 (2025), §II-A (doi:10.1109/TSP.2025.3573044, read in arXiv:2505.15318v1). The question: does PnP-ISTA with a row-normalised kernel denoiser converge for every admissible step size?

Section numbers below refer to `gaussian_pnp_stability.pdf`.

* **Proved in writing**
  * Theorem 1(a), n ≤ 3. Irreducible unit-diagonal kernels with entries at most one, and B ⪰ 0 with Be = be, give Schur-stable PnP-ISTA and PnP-ADMM for all 0 < t < 2/ρ(B). The proof reduces to Banerjee–Chaudhury Theorem 1.3 via the positivity identity (6)–(7), and the paper also gives a self-contained argument.
  * Theorem 2: a gap-dependent bound on the safe step.
  * Theorem 5: the three-state reentrant family without unit diagonals.
  * Theorem 6 and Corollary 1: symmetric denoisers. This extends Athalye, Chaudhury and Kumar (2023).
  * Lemma 1.
* **Proved with exact rational certificates** (certificate data in Appendix C)
  * Theorem 1(b): an explicit four-state Gaussian family that is unstable for each prescribed 0 < t ≤ 1/2, with a fixed positive forward operator.
  * Theorem 4: an exact two-parameter classification on [0, 10⁻⁴] × (0, 2), with ε_crit ≈ 1.6947×10⁻⁵ and t* ≈ 0.5313.
  * Theorem 3: a twelve-state positive circulant (convolution) family, unstable exactly for 0 < t < t_c ≈ 1.9157×10⁻⁴ along ε = t².
* **Numerical illustration only**: the Section 12 growth table and Figure 1.
* **Open**: everything in Section 14. This includes large steps in dimension four, whether √κ_π is needed, the case Be ≤ ρ(B)e, minimal circulant dimension and natural-image patch realisations.

**Practical size.** The divergence is extremely slow, at most about 1.7×10⁻⁵ per iteration in the four-state family and about 3.6×10⁻⁸ in the circulant family. The affinities involved are tiny, between 10⁻⁸ and 10⁻⁶. The results are structural and do not demonstrate divergence in an imaging pipeline.

**What is not claimed.**
* The first disproof of Banerjee–Chaudhury Conjecture 1.1. That counterexample is S. Patel's MathDB post. MathDB labels it unverified; a referee of this paper checked it.
* Any imaging experiment.
* Minimality of twelve for circulants.
* Historical priority. The literature search was bounded.
