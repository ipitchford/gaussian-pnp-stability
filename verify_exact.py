"""Exact four-dimensional Gaussian kernel instability calculation.

Exact checks for Theorem 1(b) and eq. (33) of gaussian_pnp_stability.md.
All determinant coefficients use exact rational arithmetic; no third-party
Python packages are needed. Run with python3 verify_exact.py.
The script prints the complete reduced determinant certificate and independently
checks values using the full 4-by-4 determinant.
"""

from fractions import Fraction as F
from itertools import permutations
from math import log, sqrt, exp


def matrix(rows):
    return [[F(x) for x in row] for row in rows]


def identity(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def add(a, b):
    return [[x + y for x, y in zip(r, s)] for r, s in zip(a, b)]


def scale(c, a):
    return [[c * x for x in row] for row in a]


def multiply(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def diagonal(values):
    return [[values[i] if i == j else F(0)
             for j in range(len(values))] for i in range(len(values))]


def poly_add(p, q):
    return [(p[i] if i < len(p) else F(0))
            + (q[i] if i < len(q) else F(0))
            for i in range(max(len(p), len(q)))]


def poly_multiply(p, q):
    result = [F(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            result[i + j] += x * y
    return result


def poly_det(a):
    result = [F(0)]
    n = len(a)
    for permutation in permutations(range(n)):
        sign = (-1) ** sum(permutation[i] > permutation[j]
                           for i in range(n) for j in range(i + 1, n))
        term = [F(sign)]
        for i, j in enumerate(permutation):
            term = poly_multiply(term, a[i][j])
        result = poly_add(result, term)
    return result


def determinant(a):
    return poly_det([[[x] for x in row] for row in a])[0]


def evaluate(p, t):
    return sum(x * t ** i for i, x in enumerate(p))


def aggregate(a):
    return [[a[i][0] + a[i][1], a[i][2], a[i][3]] for i in [0, 2, 3]]


I = identity(4)
J = matrix([[1] * 4] * 4)
v = [-8, -8, 17, -1]
V = matrix([[x * y for y in v] for x in v])
A = scale(F(1, 1000), add(add(scale(20, I), scale(245, J)), V))
B = multiply(A, A)
assert all(x > 0 for row in A for x in row)
assert all(x > 0 for row in B for x in row)
assert [sum(row) for row in A] == [F(1)] * 4
assert [sum(row) for row in B] == [F(1)] * 4
assert min(x for row in A for x in row) == F(109, 1000)
assert B == add(add(scale(F(1, 2500), I), scale(F(2499, 10000), J)),
                scale(F(229, 500000), V))
K0 = matrix([[1, F(9, 10), F(1, 10), 0],
             [F(9, 10), 1, F(1, 10), 0],
             [F(1, 10), F(1, 10), 1, 0],
             [0, 0, 0, 1]])
H = matrix([[0, F(3, 2), F(3, 10), F(1, 5)],
            [F(3, 2), 0, F(3, 10), F(1, 5)],
            [F(3, 10), F(3, 10), 0, F(3, 5)],
            [F(1, 5), F(1, 5), F(3, 5), 0]])
d0 = [sum(row) for row in K0]
assert [sum(row) for row in H] == d0
eta = F(1, 100000)

Kb0, Hb, Bb = map(aggregate, [K0, H, B])
Db0 = diagonal([sum(row) for row in Kb0])
Lb0 = add(Db0, scale(-1, Kb0))
Lb1 = add(Db0, scale(-1, Hb))
Kb0B = multiply(Kb0, Bb)
HbB = multiply(Hb, Bb)
Mpoly = [[[Lb0[i][j], Kb0B[i][j] + eta * Lb1[i][j],
           eta * HbB[i][j]] for j in range(3)] for i in range(3)]
p = poly_det(Mpoly)
f = [F(-1973, 25000000), F(74272809, 500000000000),
     F(-236896251, 625000000000000000),
     F(-4748139, 1562500000000000000000)]
assert p == [F(0), F(0)] + poly_multiply([F(1), eta], f)
assert f[1] < F(3, 20000)
assert f[2] < 0 and f[3] < 0
assert f[0] + F(1, 2) * F(3, 20000) == F(-49, 12500000)
print('PASS: universal negativity certificate for every 0 < t <= 1/2')
print('reduced determinant coefficients:', list(map(str, p)))
print('cubic after factoring t^2*(1+t/100000):', list(map(str, f)))
print('fixed stationary distribution:', [str(x / sum(d0)) for x in d0])

for t in [F(1, 1000000), F(1, 10000), F(1, 100), F(1, 10), F(1, 2)]:
    epsilon = eta * t
    K = add(K0, scale(epsilon, H))
    D = diagonal([sum(row) for row in K])
    M = add(add(D, scale(-1, K)), scale(t, multiply(K, B)))
    W = [[K[i][j] / D[i][i] for j in range(4)] for i in range(4)]
    P = multiply(W, add(I, scale(-t, B)))
    antisymmetric_M = (2 * (1 + epsilon)
                      - (F(1, 10) - F(3, 2) * epsilon) * (1 - t / 2500))
    reduced_value = evaluate(p, t)
    assert determinant(M) == antisymmetric_M * reduced_value
    assert determinant(add(I, scale(-1, P))) < 0
    reverse_numerator = add(add(I, scale(-1, W)), scale(t, multiply(B, W)))
    assert determinant(reverse_numerator) == determinant(add(I, scale(-1, P)))
    assert determinant(add(I, scale(t, B))) > 0
    assert reduced_value < 0
    a = float(K[0][1])
    b = float(K[0][2])
    u, w = -log(a), -log(b)
    h = sqrt(w - u / 4)
    y = (w - u / 2 + log(3)) / (2 * h)
    z = sqrt(-log(float(epsilon / 5)) - u / 4 - y * y)
    points = [(-sqrt(u) / 2, 0, 0), (sqrt(u) / 2, 0, 0),
              (0, h, 0), (0, y, z)]
    max_error = max(abs(exp(-sum((points[i][k] - points[j][k]) ** 2
                                 for k in range(3))) - float(K[i][j]))
                    for i in range(4) for j in range(4))
    assert max_error < 1e-14
    print('verified t =', t, '; Gaussian entry error =', max_error)


# An independent two-parameter calculation checks the first-loss asymptotics.
# Derive the quotient matrices from the original four-dimensional data rather
# than entering the quotient matrices from the manuscript by hand.
def difference_quotient(a):
    return [[a[i][j] - a[2][j] for j in range(2)] for i in range(2)]


def sparse_add(a, b):
    result = dict(a)
    for key, value in b.items():
        result[key] = result.get(key, F(0)) + value
    return {key: value for key, value in result.items() if value}


def sparse_multiply(a, b):
    result = {}
    for (ei, ti), x in a.items():
        for (ej, tj), y in b.items():
            key = (ei + ej, ti + tj)
            result[key] = result.get(key, F(0)) + x * y
    return {key: value for key, value in result.items() if value}


degree_bar = [sum(row) for row in Kb0]
Q0 = difference_quotient([[x / degree_bar[i] for x in row]
                          for i, row in enumerate(Kb0)])
Q1 = difference_quotient([[x / degree_bar[i] for x in row]
                          for i, row in enumerate(Hb)])
C = difference_quotient(Bb)
Q0C, Q1C = multiply(Q0, C), multiply(Q1, C)
I2 = identity(2)
# Dictionary indices are (power of epsilon, power of t).
two_parameter = [[{(0, 0): I2[i][j] - Q0[i][j],
                   (1, 0): I2[i][j] - Q1[i][j],
                   (0, 1): Q0C[i][j],
                   (1, 1): Q1C[i][j]}
                  for j in range(2)] for i in range(2)]
term1 = sparse_multiply(two_parameter[0][0], two_parameter[1][1])
term2 = sparse_multiply(two_parameter[0][1], two_parameter[1][0])
delta = sparse_add(term1, {key: -value for key, value in term2.items()})
expected_delta = {
    (1, 0): F(31, 120), (2, 0): F(217, 200),
    (0, 1): F(-133, 3750000), (1, 1): F(1782779, 10000000),
    (2, 1): F(969583, 50000000), (0, 2): F(751389, 12500000000),
    (1, 2): F(-15987, 1000000000),
    (2, 2): F(-1582713, 125000000000),
}
assert delta == expected_delta
assert sum(Q0[i][i] for i in range(2)) == 1 + F(47, 60)
assert determinant(Q0) == F(47, 60)
fast_gap = 1 - F(47, 60)
gap_derivative = delta[(1, 0)] / fast_gap
step_derivative = -delta[(0, 1)] / fast_gap
threshold_derivative = gap_derivative / step_derivative
assert gap_derivative == F(31, 26)
assert step_derivative == F(133, 812500)
assert threshold_derivative == F(968750, 133)
assert threshold_derivative / gap_derivative == F(812500, 133)
assert step_derivative - eta * gap_derivative == F(1973, 13000000)
print('PASS: all eight bivariate determinant coefficients')
print('gap derivative:', gap_derivative)
print('first-loss derivative:', threshold_derivative)
print('threshold / gap limit:', threshold_derivative / gap_derivative)
print('All exact checks passed. Gaussian floating-point checks are supplementary.')
