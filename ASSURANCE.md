# Assurance

This file sets out what has been checked and what has not. Each row is a separate kind of assurance, and no row implies another.

| Dimension | State | Evidence |
|---|---|---|
| Availability | Public GitHub release and Zenodo DOI 10.5281/zenodo.23690125 | Release assets, [SHA256SUMS](SHA256SUMS) |
| Written proofs | Complete proofs in the manuscript | `gaussian_pnp_stability.pdf`, Sections 5–11. Theorem 1(a) uses Banerjee–Chaudhury Theorem 1.3 as published. A self-contained argument is also given. |
| Exact certificates | Exact rational arithmetic | `verify_exact.py`, `verify_circulant.py` and `verify_reentrant.py` (standard library); `certify_revision.py` (sympy) |
| Certificates checkable from the paper | Printed in Appendix C and in hashed coefficient files | Appendix C shows the polynomials, Bernstein sign summaries, Sturm counts and isolating intervals for ε_crit, t* and t_c. The data files are `certificates/*.json`. |
| Stand-alone certificate check | Passed | `check_certificates.py` uses the standard library only. It reads the coefficient files and recomputes every Bernstein array and Sturm chain. Output: `certificate_check_output.txt`. |
| Numerical results | Labelled as such | The Section 12 table (`practical_numbers.py`, 50-digit floating point) and Figure 1 |
| Internal replay | Passed | `replay_all.py` reruns every verifier from a temporary copy and compares the records. See [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md). |
| Semantic negative controls | Passed | `negative_controls.py` applies eight corruptions to the certificate files: a Bernstein sign, a side-condition polynomial, a Δ coefficient, the ε_crit interval, the t* interval, a t* coefficient, the t_c interval and the sign of p(0). The checker must reject each one. |
| Optimised mode | Not applicable by design | The exact checks use `assert`, so they must not run under `python -O`. `check_certificates.py` and the negative controls do not rely on `assert`. |
| Implementation diversity | Partial | For a targeted proof audit, separate code was written that rebuilt the Theorem 3 and 4 polynomials from the matrix definitions. It used its own Bernstein subdivision and root isolation, and reproduced the t* and t_c isolating intervals. This is diversity within one producer workflow, not independent reimplementation. |
| Internal AI-assisted review | Done | A five-perspective review, AI-assisted correctness audits by a different model of the draft and of the revision, a targeted proof audit, a verification re-review and final editor fixes. See [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). This is not peer review. |
| Independent reproduction by an unaffiliated party | Not established | |
| Specialist (human) review | Not established | |
| Editorial peer review | Not established | |
| Formal verification | Not attempted | No proof assistant was used. |
| Source fidelity | Checked | Sinha et al.'s open-question sentence was quoted verbatim from arXiv:2505.15318v1, §II-A. The hypotheses of Banerjee–Chaudhury Theorems 1.2 and 1.3 and the Athalye et al. (2023) results were read in the arXiv full texts. |
| Novelty | Bounded search only | [RESEARCH_GATES.json](RESEARCH_GATES.json), `novelty_report.md` |
| Citations | Every DOI resolved and its title matched | [CITATION_AUDIT.md](CITATION_AUDIT.md) |

## Scope and limits

* The four-dimensional instability is certified only for 0 < t < t* ≈ 0.5313, with the kernel chosen for each step. Each fixed primitive kernel is stable for small enough steps (Banerjee–Chaudhury Theorem 1.1).
* The convolution instability is certified only for 0 < t < 1.92×10⁻⁴.
* The growth is at most 1.7×10⁻⁵ per iteration. The features are abstract points, not image patches.
* A byte-identical replay shows that the shipped records are what the shipped programs produce. It does not by itself show that the programs implement the manuscript's definitions. That rests on the written proofs, the printed certificate data, the stand-alone checker and the negative controls.
