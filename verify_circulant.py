"""Exact verification of the twelve-state positive circulant example.

Reproduces the constants of Theorem 3 of gaussian_pnp_stability.md and an
exact witness OFF the curve eps=t^2 (t=1/100, eps=1e-9). The certified range
and the witness ON the curve are in certify_revision.py.  Only the Python standard
library is required.  The finite determinant signs are exact; the separate
Gaussian coordinate check uses floating-point arithmetic as a sanity check.
The arbitrary-small-step assertion is proved analytically in the manuscript.
"""

from fractions import Fraction as F
from math import exp, log, sqrt
from pathlib import Path
import json


def identity(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def add(a, b):
    return [[x + y for x, y in zip(r, s)] for r, s in zip(a, b)]


def scale(c, a):
    return [[c * x for x in row] for row in a]


def multiply(a, b):
    return [[sum(x * y for x, y in zip(row, column))
             for column in zip(*b)] for row in a]


def determinant(a):
    """Exact Gaussian elimination; row pivoting prevents zero divisions."""
    a = [row[:] for row in a]
    result = F(1)
    n = len(a)
    for i in range(n):
        pivot = next((j for j in range(i, n) if a[j][i]), None)
        if pivot is None:
            return F(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            result = -result
        value = a[i][i]
        result *= value
        for j in range(i + 1, n):
            coefficient = a[j][i] / value
            for k in range(i, n):
                a[j][k] -= coefficient * a[i][k]
    return result


n = 12
group_g = [0, 1, 2, 4, 6]
group_h = [3, 5, 7, 8, 9, 10, 11]
high_g = {0, 2, 4, 6}
high_h = {5, 8, 9, 10}
epsilon, t = F(1, 10**9), F(1, 100)
I = identity(n)
J = [[F(1)] * n for _ in range(n)]
E = scale(F(1, n), J)
Z = [[F([1, 0, -1, 0][(i - j) % 4], 6)
      for j in range(n)] for i in range(n)]
A = add(add(scale(F(11, 12), E), scale(F(1, 12), I)),
        scale(F(5, 12), Z))
B = multiply(A, A)
assert multiply(Z, Z) == Z
assert all(sum(row) == 0 for row in Z)
assert B == add(add(scale(F(143, 144), E), scale(F(1, 144), I)),
                scale(F(35, 144), Z))
assert all(x > 0 for row in A for x in row)
assert min(x for row in A for x in row) == F(1, 144)
assert all(sum(row) == 1 for row in A)
assert all(A[i][j] == A[0][(j - i) % n]
           for i in range(n) for j in range(n))

K0 = [[F(0)] * n for _ in range(n)]
for group, high in [(group_g, high_g), (group_h, high_h)]:
    for i in group:
        for j in group:
            K0[i][j] = (F(1) if i == j else F(99, 100)
                        if i in high and j in high else F(1, 100))
K = add(scale(1 - epsilon, K0), scale(epsilon, J))
degrees = [sum(row) for row in K]
W = [[x / degrees[i] for x in row] for i, row in enumerate(K)]
pi = [degree / sum(degrees) for degree in degrees]
assert all(x > 0 for row in K for x in row)
assert all(K[i][i] == 1 for i in range(n))
assert all(W[i][j] * pi[i] == W[j][i] * pi[j]
           for i in range(n) for j in range(n))

P = multiply(W, add(I, scale(-t, B)))
det_p = determinant(add(I, scale(-1, P)))
denominator_r = determinant(add(I, scale(t, B)))
assert denominator_r == (1 + t) * (1 + t / 4)**2 * (1 + t / 144)**9
reverse_numerator = add(add(I, scale(-1, W)), scale(t, multiply(B, W)))
assert determinant(reverse_numerator) == det_p
det_r = det_p / denominator_r
assert det_p < 0 and det_r < 0

# Compute stationary block weights directly from K0, then compress B and Z.
degrees0 = [sum(row) for row in K0]
weight_g = [degrees0[i] / sum(degrees0[j] for j in group_g) for i in group_g]
weight_h = [degrees0[i] / sum(degrees0[j] for j in group_h) for i in group_h]
assert weight_g == [F(199, 848) if i in high_g else F(13, 212)
                    for i in group_g]
assert weight_h == [F(200, 959) if i in high_h else F(53, 959)
                    for i in group_h]
c_z = (sum(weight_g[k] * sum(Z[i][j] for j in group_g)
           for k, i in enumerate(group_g))
       - sum(weight_h[k] * sum(Z[i][j] for j in group_g)
             for k, i in enumerate(group_h)))
c_b = F(1, 144) + F(35, 144) * c_z
assert c_z == F(-12875, 406616)
assert c_b == F(-6287, 8364672)

# Hierarchical Gaussian features: 12 private + 2 shared + 1 block direction.
a = F(99, 100) + epsilon / 100
b = F(1, 100) + F(99, 100) * epsilon
h, ell, distance = -log(float(a)), -log(float(b)), -log(float(epsilon))
points = []
for i in range(n):
    point = [0.0] * 15
    if i in high_g or i in high_h:
        point[i] = sqrt(h / 2)
        point[12 if i in high_g else 13] = sqrt((ell - h) / 2)
    else:
        point[i] = sqrt(ell / 2)
    point[14] = (1 if i in group_g else -1) * sqrt(distance - ell) / 2
    points.append(point)
maximum_error = max(
    abs(exp(-sum((x - y)**2 for x, y in zip(points[i], points[j])))
        - float(K[i][j])) for i in range(n) for j in range(n)
)
assert maximum_error < 1e-14

certificate = {
    'n': n,
    'epsilon': epsilon,
    't': t,
    'forward_A_first_row': A[0],
    'B_first_row': B[0],
    'kernel_blocks': [group_g, group_h],
    'high_similarity_groups': [sorted(high_g), sorted(high_h)],
    'K0_high_similarity': F(99, 100),
    'K0_other_inblock_similarity': F(1, 100),
    'c_projection': c_z,
    'c_B': c_b,
    'det_I_P': det_p,
    'det_I_R': det_r,
    'det_I_plus_tB': denominator_r,
    'all_exact_checks_passed': True,
}
destination = Path(__file__).with_name('circulant_certificate.json')
destination.write_text(json.dumps(certificate, default=str, indent=2) + '\n')
print('PASS: positive symmetric circulant forward matrix, condition(A)=12')
print('PASS: positive unit-diagonal reversible kernel')
print('PASS: exact det(I-P)<0 and det(I-R)<0 at t=1/100, epsilon=10^-9')
print('compressed projection eigenvalue:', c_z)
print('compressed Hessian eigenvalue:', c_b)
print('Gaussian coordinate error (supplementary):', maximum_error)
print('Exact finite certificate written to', destination.name)
