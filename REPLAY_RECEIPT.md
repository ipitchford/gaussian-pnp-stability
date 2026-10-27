# Replay receipt

This is a producer-side replay, run on 27 October 2026 on macOS (Apple silicon), using Python 3.11.14 via `uv run python3 -I -B`, from the assembled release tree. It is not an independent reproduction.

| Command | Result |
|---|---|
| `python3 -I -B replay_all.py` (sympy 1.13.3) | PASS. The run copies the package to a temporary directory and runs eight scripts there: `verify_exact`, `verify_circulant`, `verify_reentrant`, `certify_revision`, `emit_certificates`, `check_certificates`, `practical_numbers` and `make_fingerprints`. All exit 0. The three certificate files, `circulant_certificate.json` and `exact_fingerprints.json` are byte-identical to the shipped copies. `revision_certificate.json` and `reentrant_verification.json` are identical apart from timing and version fields. The checker's stdout equals `certificate_check_output.txt`. Wall time is about 40 s. |
| `python3 -I -B negative_controls.py` | PASS. The unmodified certificates are accepted, and all 8 corruptions are rejected. |
| `python3 -I -B check_certificates.py` | PASS. All Bernstein arrays are positive, and the Sturm counts and isolating intervals for ε_crit, t* and t_c hold. Output: `certificate_check_output.txt`. |
| `python3 build_pdfs.py` | 28 pages. Every page was inspected, and a PDF TeX preflight found no raw TeX. |

**Earlier replays.** Before packaging, every script, including `make_figure.py`, was rerun twice from fresh copies (`verification_output.txt`). The figure PNG was byte-identical on this machine. The targeted proof audit independently rebuilt the Theorem 3 and 4 polynomials from the matrix definitions.

**Continuous integration.** GitHub Actions (`.github/workflows/replay.yml`) repeats the following on Ubuntu with Python 3.11: the checksum check, the checker, the negative controls and the full replay.
