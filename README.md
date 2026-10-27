# Gaussian-kernel plug-and-play iterations can be unstable

**A sharp dimension threshold for irreducible unit-diagonal kernels, the role of the mixing gap, and why unit diagonals matter.**

**Release 1.0.0-candidate (27 October 2026). Unrefereed research candidate.**

- DOI: [10.5281/zenodo.23690125](https://doi.org/10.5281/zenodo.23690125) (this version)
- Concept DOI: [10.5281/zenodo.23690124](https://doi.org/10.5281/zenodo.23690124) (all versions)
- Repository: <https://github.com/ipitchford/gaussian-pnp-stability> (tag `v1.0.0-candidate`)
- Evidence Press page: <https://evidencepress.org/releases/gaussian-pnp-stability/>

The package was checked by internal AI-assisted review only. That is not peer review, independent reproduction or formal verification. See [STATUS.md](STATUS.md) and [ASSURANCE.md](ASSURANCE.md).

Author: Anonymous. Text and data are CC0-1.0; code is MIT (see `LICENSE`).

## The question

Sinha, Kumar, Athalye and Chaudhury (IEEE Trans. Signal Process. 73, 2025, §II-A) observed that PnP-ISTA with row-normalised kernel denoisers, such as Gaussian non-local means, seems to converge in practice. They said they had "not yet found a formal proof or a counterexample". This package answers the question for PnP-ISTA and PnP-ADMM.

The manuscript is `gaussian_pnp_stability.pdf` (28 pages). Its source is `gaussian_pnp_stability.md`.

## Main results

These hold for irreducible unit-diagonal kernels with entries at most one, and a data Hessian B ⪰ 0 with Be = be:

- **Dimensions 1–3 are stable** for every step 0 < t < 2/ρ(B) (Theorem 1(a)). The proof reduces to Banerjee–Chaudhury Theorem 1.3 via a new positivity identity.
- **Dimension 4 can be unstable.** Explicit Gaussian kernels and a fixed positive forward operator make both iterations diverge at any prescribed step 0 < t ≤ 1/2 (Theorem 1(b)). Theorem 4 classifies the family exactly. Its boundary is t* ≈ 0.5313.
- **Convolution operators are not immune.** A twelve-state positive circulant example fails exactly for 0 < t < t_c ≈ 1.9157×10⁻⁴ (Theorem 3).
- **The cause is a small mixing gap g.** The step t < g/√κ_π is always safe (Theorem 2), and linear scaling in g cannot be improved.
- **Unit diagonals matter.** Without them, failure occurs already in dimension 3 (Theorem 5).
- **Symmetric denoisers** are stable exactly when ker(I−W) ∩ ker B = {0} (Theorem 6, extending Athalye, Chaudhury and Kumar, 2023). Symmetrically normalised irreducible kernels always satisfy this when Be = be (Corollary 1).

**How big is the effect? Tiny.** The divergence is at most about 1.7×10⁻⁵ per iteration, so tenfold growth takes 10⁵–10⁶ iterations, against the 20–500 of a typical run. The affinities involved are 10⁻⁸ to 10⁻⁶. The results show that no uniform step-size guarantee exists. They do not show divergence in a real imaging pipeline.

**Practical advice.** Scale the step with the denoiser's mixing gap, or symmetrise the denoiser (Section 12).

## Reproduce

The code was tested with Python 3.11 via uv. Run the scripts in a writable copy, because they write their outputs beside themselves.

```text
sha256sum -c SHA256SUMS                                   # or: shasum -a 256 -c SHA256SUMS
uv run python3 -I -B check_certificates.py                # Appendix C, standard library, seconds
uv run python3 -I -B negative_controls.py                 # 8 corrupted certificates must be rejected
uv run --with sympy==1.13.3 python3 -I -B replay_all.py   # every verifier from a temporary copy, about 1 min
```

Individual scripts:

```text
uv run python3 -I verify_exact.py
uv run python3 -I verify_circulant.py
uv run python3 -I verify_reentrant.py
uv run --with sympy==1.13.3 python3 -I certify_revision.py
uv run --with sympy==1.13.3 python3 -I emit_certificates.py
uv run python3 -I check_certificates.py
uv run --with sympy==1.13.3 python3 -I practical_numbers.py
uv run --with numpy==2.1.3 --with matplotlib==3.9.2 python3 -I make_figure.py
uv run python3 -I make_fingerprints.py
python3 build_pdfs.py    # needs pandoc and pdfLaTeX
```

Figures are reproducible visually, but not byte for byte across matplotlib builds.

## Files

**Manuscript**
- `gaussian_pnp_stability.md` / `.pdf`: the manuscript. It includes the merged unit-diagonal supplement (Section 10) and the certificate appendix (Appendix C).

**Verifiers**
- `verify_exact.py`: the certificate for Theorem 1(b) (standard library).
- `verify_circulant.py`, `circulant_certificate.json`: the constants of Theorem 3 and an exact off-curve witness.
- `verify_reentrant.py`, `reentrant_verification.json`, `reentrant_verification_output.txt`: Theorem 5.
- `certify_revision.py`, `revision_certificate.json`: Theorem 4, the certified range of Theorem 3 with an on-curve witness, and the conditioning bound 19.22 (needs sympy).

**Certificates for Theorems 3 and 4 (Appendix C)**
- `emit_certificates.py` (sympy) writes the exact coefficient files in `certificates/*.json`.
- `check_certificates.py` (standard library only) re-verifies every Bernstein sign, Sturm count and isolating interval from those files. Its output is in `certificate_check_output.txt`.

**Numbers and figure**
- `practical_numbers.py`: the growth rates of Section 12 (50-digit floating point; illustration only).
- `make_figure.py`, `gaussian_stability_diagram.png/.pdf`: Figure 1.

**Replay**
- `replay_all.py`, `negative_controls.py`: the replay driver and the semantic negative controls.
- `verification_output.txt`: captured output of every script from a clean copy.

**Pre-revision records and build**
- `novelty_report.md`, `proof_audit.md`, `proof_audit_verification.txt`, `independent_verification.md`, `exact_fingerprints.json`, `make_fingerprints.py`: pre-revision records (20 October 2026), superseded where they differ from the manuscript.
- `target.yaml`: the claim record, updated to the final version.
- `build_pdfs.py`, `pdf_header.tex`: the PDF build.

## Release records

- [STATUS.md](STATUS.md)
- [ASSURANCE.md](ASSURANCE.md)
- [CLAIMS.json](CLAIMS.json)
- [AI_INDEX.md](AI_INDEX.md)
- [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md)
- [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md)
- [ENVIRONMENT.txt](ENVIRONMENT.txt)
- [SOURCES.md](SOURCES.md)
- [CITATION_AUDIT.md](CITATION_AUDIT.md)
- [RESEARCH_GATES.json](RESEARCH_GATES.json)
- [PROVENANCE.md](PROVENANCE.md)
- [LICENSES.md](LICENSES.md)
- [PUBLIC_DOMAIN.md](PUBLIC_DOMAIN.md)

`SHA256SUMS` lists every file.

## Citation

Anonymous (2026). *Gaussian-kernel plug-and-play iterations can be unstable: a sharp dimension threshold for irreducible unit-diagonal kernels, the role of the mixing gap, and why unit diagonals matter* (Version 1.0.0-candidate) [Unrefereed research candidate]. Evidence Press. https://doi.org/10.5281/zenodo.23690125

The scholarly creator is Anonymous. Ian Pitchford maintains and publishes the repository; that role is not authorship.
