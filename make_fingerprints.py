#!/usr/bin/env python3
"""Exact finite-object fingerprints for the bounded novelty audit.

This records object identity; it does not establish historical novelty.
Only the Python standard library is used.
"""
import contextlib
from fractions import Fraction as F
import hashlib
import io
from itertools import permutations
import json
from pathlib import Path
import runpy

HERE = Path(__file__).resolve().parent
with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(HERE/"verify_exact.py"))

def serial(matrix):
    return [[str(x) for x in row] for row in matrix]

def digest(matrix):
    return hashlib.sha256(json.dumps(serial(matrix),separators=(",",":")).encode()).hexdigest()

def fingerprint(matrix):
    n = len(matrix)
    transformed = []
    for order in permutations(range(n)):
        transformed.append(json.dumps(serial([[matrix[i][j] for j in order] for i in order]),separators=(",",":")))
    char = ns["poly_det"]([[[(-matrix[i][j]),F(i==j)] for j in range(n)] for i in range(n)])
    return {
        "matrix": serial(matrix), "dimension": n,
        "rank": n if ns["determinant"](matrix) else "not inferred",
        "trace": str(sum(matrix[i][i] for i in range(n))),
        "determinant": str(ns["determinant"](matrix)),
        "characteristic_coefficients_ascending": list(map(str,char)),
        "row_sums": [str(sum(row)) for row in matrix],
        "column_sums": [str(sum(row)) for row in zip(*matrix)],
        "sorted_entry_multiset": list(map(str,sorted(x for row in matrix for x in row))),
        "strictly_positive": all(x>0 for row in matrix for x in row),
        "symmetric": all(matrix[i][j]==matrix[j][i] for i in range(n) for j in range(n)),
        "sha256_exact_rows": digest(matrix),
        "sha256_transpose": digest([list(row) for row in zip(*matrix)]),
        "sha256_index_reversal": digest([list(reversed(row)) for row in reversed(matrix)]),
        "sha256_canonical_simultaneous_permutation": hashlib.sha256(min(transformed).encode()).hexdigest(),
        "canonicalisation": "lexicographic minimum of exact Fraction strings over all n! simultaneous row/column permutations"
    }

instances = []
for t in [F(1,100),F(1,4),F(1,2)]:
    eps = t/100000
    K = ns["add"](ns["K0"],ns["scale"](eps,ns["H"]))
    W = [[x/sum(row) for x in row] for row in K]
    instances.append({"t":str(t),"epsilon":str(eps),"K":fingerprint(K),"W":fingerprint(W)})
result = {
    "date":"2026-10-20",
    "contribution_unit":"structured stability theorem, not historical priority",
    "A":fingerprint(ns["A"]), "B":fingerprint(ns["B"]),
    "A_spectrum_with_multiplicity":["1","219/500","1/50","1/50"],
    "B_spectrum_with_multiplicity":["1","47961/250000","1/2500","1/2500"],
    "fixed_stationary_distribution":["10/31","10/31","6/31","5/31"],
    "instances":instances,
    "instability_cubic_coefficients":[str(x) for x in ns["f"]],
    "first_loss_derivative":"968750/133",
    "first_loss_over_gap_limit":"812500/133",
    "notes":["No generating-function or enumerative sequence identification is claimed.",
             "These fingerprints permit exact equality tests; absence of a search match does not establish novelty."]
}
(HERE/"exact_fingerprints.json").write_text(json.dumps(result,indent=2)+"\n")
print("Saved exact_fingerprints.json: A, B and three rational Gaussian-kernel instances.")
