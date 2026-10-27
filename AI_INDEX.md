# Agent navigation

## Identity

- Release slug: `gaussian-pnp-stability`. Version `1.0.0-candidate`, tag `v1.0.0-candidate`.
- DOI 10.5281/zenodo.23690125; concept DOI 10.5281/zenodo.23690124.
- Repository: https://github.com/ipitchford/gaussian-pnp-stability
- Release page: https://evidencepress.org/releases/gaussian-pnp-stability/
- Scholarly creator: Anonymous. Status: unrefereed candidate ([STATUS.md](STATUS.md)).
- Manuscript: [PDF](gaussian_pnp_stability.pdf) (28 pages) and [Markdown source](gaussian_pnp_stability.md). Figure: [gaussian_stability_diagram.png](gaussian_stability_diagram.png) (illustration only).
- Overview and claims: [README.md](README.md), [CLAIMS.json](CLAIMS.json), [ASSURANCE.md](ASSURANCE.md), [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md).
- Sources and provenance: [PROVENANCE.md](PROVENANCE.md), [SOURCES.md](SOURCES.md), [CITATION_AUDIT.md](CITATION_AUDIT.md), [RESEARCH_GATES.json](RESEARCH_GATES.json), [novelty_report.md](novelty_report.md) (pre-revision).
- Citation and licences: [CITATION.cff](CITATION.cff), [.zenodo.json](.zenodo.json), [LICENSES.md](LICENSES.md), [PUBLIC_DOMAIN.md](PUBLIC_DOMAIN.md).

## Scientific map

**Setting.** The two iteration matrices are P(t) = W(I − tB) for PnP-ISTA and R(t) = I − W + (I + tB)⁻¹(2W − I) for PnP-ADMM, where:
- W = diag(Ke)⁻¹K is the row-normalised kernel denoiser;
- K is a symmetric nonnegative affinity matrix, Gaussian when K_ij = exp(−|x_i − x_j|²);
- B = AᵀA ⪰ 0.

"Stable" means Schur stable (ρ < 1). The standing class is K irreducible with K_ii = 1 and K_ij ≤ 1, and Be = be with b > 0. The target is the open question in Sinha et al. (2025), §II-A.

**Results.**
- n ≤ 3: stable for all 0 < t < 2/ρ(B). Theorem 1(a), **proved**, via Banerjee–Chaudhury Theorem 1.3 and the identity (6)–(7).
- n = 4: an explicit Gaussian family is unstable for each 0 < t ≤ 1/2. Theorem 1(b), **exact certificate**, [verify_exact.py](verify_exact.py).
- An exact classification of that family, with ε_crit ≈ 1.6947×10⁻⁵ and t* ≈ 0.5313. Theorem 4, **exact certificate**: [certificates/thm4_classification.json](certificates/thm4_classification.json), [certificates/thm4_tstar.json](certificates/thm4_tstar.json).
- A twelve-state positive circulant family is unstable exactly for 0 < t < t_c ≈ 1.9157×10⁻⁴ along ε = t². Theorem 3, **exact certificate**: [certificates/thm3_circulant.json](certificates/thm3_circulant.json).
- The step t < g/√κ_π is safe. Theorem 2, **proved**. The order g is sharp in the family, but the constant is not.
- Without unit diagonals, n = 3 fails. Theorem 5, **proved**, [verify_reentrant.py](verify_reentrant.py).
- Symmetric W is stable iff ker(I − W) ∩ ker B = {0}. Theorem 6, **proved**, extending Athalye, Chaudhury and Kumar (2023).

**Do not read:**
- the instability as practically large: it is at most 1.7×10⁻⁵ per iteration, and the affinities are 10⁻⁸ to 10⁻⁶;
- any result as an imaging experiment;
- this as the first disproof of Banerjee–Chaudhury Conjecture 1.1 (that was Patel's MathDB post);
- Theorem 1(b) as a single fixed kernel unstable at all small steps;
- internal AI-assisted review as peer review.

## Code and replay

[check_certificates.py](check_certificates.py) uses the standard library and needs seconds. [negative_controls.py](negative_controls.py) runs 8 controls. [replay_all.py](replay_all.py) needs sympy 1.13.3. [emit_certificates.py](emit_certificates.py) builds the coefficient files from the matrix definitions.

Results are in [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), [verification_output.txt](verification_output.txt) and [certificate_check_output.txt](certificate_check_output.txt). CI: [.github/workflows/replay.yml](.github/workflows/replay.yml).

## Useful next research question

Is any irreducible unit-diagonal Gaussian kernel on four points unstable at a step t ∈ (t*, 2) (§14.1)? Can the obstruction be realised with natural-image patches inside a search window (§14.6)?
