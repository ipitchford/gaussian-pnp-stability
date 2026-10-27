# Review summary (internal, AI-assisted)

This file summarises the substantive objections raised in internal AI-assisted review and how each was resolved. It is not peer review. The original review records and the point-by-point response are kept privately.

## Review rounds

1. A five-perspective review, from editor, methodology, domain, perspective and devil's-advocate viewpoints.
2. An AI-assisted correctness audit of the draft by a different model.
3. Revision, including merging the unit-diagonal supplement in as Section 10.
4. A verification re-review of the revision.
5. A targeted proof audit, using independent exact rebuilds of the Theorem 3 and 4 polynomials with its own Bernstein and root-isolation code.
6. A second correctness audit of the revised text, by the same different model.
7. Final editor fixes (27 October 2026).

## Objections and resolutions

| Objection | Resolution |
|---|---|
| The sharp-dimension class must include irreducibility and Be = be. | Irreducibility is now stated throughout, with a reducible n = 2 counterexample (K = I₂). |
| The step range of Theorem 1(b) was overstated. | It is restricted to t ≤ 1/2. The exact boundary t* ≈ 0.5313 is given, with stability on (t*, 2) (Theorem 4). |
| The first loss of stability was underproved. | Theorem 4 now gives an exact two-parameter classification, using Lemma 1 (the Jury test), Bernstein sign certificates and Sturm counts. |
| The certificates for Theorems 3 and 4 could not be checked from the paper (second audit, main point). | Appendix C prints the quadratic factors, the Bernstein sign summaries, the Sturm counts and rational isolating intervals for ε_crit, t* and t_c. It also gives the SHA-256 hashes of the coefficient files. A standard-library checker and eight negative controls are shipped. |
| The small-t window of Theorem 3 was not quantified. | It is certified as 0 < t < t_c = 1.9157074941…×10⁻⁴, with an exact witness on the curve. A digit error in t_c (…940) was found and corrected in the final fixes. |
| Positivity after the threshold was asserted up to t = 2, outside the Gaussian domain (second audit). | It is now stated for t_c < t < 1. The Sturm count excludes any further root up to 2, and t ≥ 1 is flagged as an algebraic continuation. |
| "Symmetric denoisers are always stable" is false without a condition (second audit). | The condition ker(I−W) ∩ ker B = {0} is now stated everywhere, with the counterexample W = I₂, B = eeᵀ/2. |
| Credit for Banerjee–Chaudhury Theorems 1.2 and 1.3. | Theorem 1.3 is credited, with the reduction in Section 5.2. Theorem 1.2 is adapted from Athalye, Chaudhury and Kumar (2023), which is now cited. Their Lemma 1 and Theorems 1–2 already give the iff characterisation for symmetric stochastic PSD W, so the novelty of Theorem 6 is narrowed to the non-stochastic and indefinite extension and to Corollary 1. |
| Patel's prior counterexample to Conjecture 1.1. | It is described accurately and credited first. The paper says the post was checked by a referee, not that it was "verified". |
| Theorem 2: constant gap, √κ_π and the hypothesis on B. | A conservative factor of about 8,600 is stated. The hypothesis is relaxed (B need not be stochastic), and a norm-comparison step is added. |
| Practical meaning. | Section 12 gives the growth rates, the size of the affinities, step-size advice and an imaging mapping. |
| Boundary statements: the neutral step at ε_crit, the ε = 0 range, the converse of Lemma 1, the n ≥ 2 case of Corollary 1, irreducibility vocabulary (second audit). | All corrected in the final fixes. |
| Revision-history phrasing in the text; min{T₊, 2}. | Removed; simplified to T₊ < 2, now certified by Δ(ε, 2) > 0. |

## Not done

Two items were declined:
- a natural-image or non-local-means experiment;
- an exact symbolic check of one floating-point identity.

Both are disclosed as non-claims. Liu, Xing and Sun (2021) was read in abstract only.
