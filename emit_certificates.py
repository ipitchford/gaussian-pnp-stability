#!/usr/bin/env python3
"""Emit the decisive exact data behind Theorems 3 and 4 as coefficient files.

Writes, into the directory certificates/ beside this file (or --dir DIR):
  thm4_classification.json  T_P, D_P, T_R, D_R, theta; numerators and
                            denominators of the side conditions and their
                            Bernstein coefficient arrays on the unit square
                            after eps = EPS1*u, t = 2*s; the discriminant,
                            a1, a2, 4*a2 + a1; isolating interval for eps_crit.
  thm4_tstar.json           integer cubic for t* and its isolating interval.
  thm3_circulant.json       integer polynomial p(t) = N_c(t)/t^2 (cleared),
                            isolating interval for t_c (the exact on-curve witness is in revision_certificate.json).
Every number is an exact rational written as a string. check_certificates.py
re-verifies all signs and root counts from these files alone, using only the
Python standard library. Requires sympy (pinned in README).
"""
import json
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix

HERE = Path(__file__).resolve().parent
OUTDIR = HERE / 'certificates'
if '--dir' in sys.argv:
    OUTDIR = Path(sys.argv[sys.argv.index('--dir') + 1])
OUTDIR.mkdir(parents=True, exist_ok=True)

eps, t, u, s = sp.symbols('epsilon t u s')
EPS1 = sp.Rational(1, 10000)
TMAX = sp.Integer(2)


def q(x):
    return str(sp.Rational(x))


def poly2_terms(expr):
    P = sp.Poly(sp.expand(expr), eps, t)
    return [[int(i), int(j), q(c)] for (i, j), c in sorted(P.terms())]


def poly1_coeffs(P):
    """Ascending integer coefficients of a univariate Poly with rational coefficients."""
    _, P = P.clear_denoms()
    return [str(int(c)) for c in reversed(P.all_coeffs())]


def bernstein(expr):
    Q = sp.Poly(sp.expand(sp.expand(expr).subs({eps: EPS1 * u, t: TMAX * s}, simultaneous=True)), u, s)
    m, n = max(Q.degree(u), 0), max(Q.degree(s), 0)
    c = dict(zip(Q.monoms(), Q.coeffs()))
    arr = [[q(sum(sp.binomial(k, i) / sp.binomial(m, i) * sp.binomial(l, j) / sp.binomial(n, j)
                  * c.get((i, j), 0) for i in range(k + 1) for j in range(l + 1)))
            for l in range(n + 1)] for k in range(m + 1)]
    return m, n, arr


def signed_fraction(expr):
    num, den = sp.fraction(sp.factor(sp.together(expr)))
    num, den = sp.expand(num), sp.expand(den)
    if sp.Poly(den, eps, t).eval({eps: 0, t: 1}) < 0:
        num, den = -num, -den
    return num, den


# ---------------------------------------------------------------- Theorem 4
v = sp.Matrix([-8, -8, 17, -1])
A = (20 * sp.eye(4) + 245 * sp.ones(4, 4) + v * v.T) / 1000
B = A * A
K0 = sp.Matrix([[1, sp.Rational(9, 10), sp.Rational(1, 10), 0],
                [sp.Rational(9, 10), 1, sp.Rational(1, 10), 0],
                [sp.Rational(1, 10), sp.Rational(1, 10), 1, 0],
                [0, 0, 0, 1]])
H = sp.Matrix([[0, sp.Rational(3, 2), sp.Rational(3, 10), sp.Rational(1, 5)],
               [sp.Rational(3, 2), 0, sp.Rational(3, 10), sp.Rational(1, 5)],
               [sp.Rational(3, 10), sp.Rational(3, 10), 0, sp.Rational(3, 5)],
               [sp.Rational(1, 5), sp.Rational(1, 5), sp.Rational(3, 5), 0]])
FF = sp.QQ.frac_field(eps, t)
K = K0 + eps * H
W = sp.diag(*[sum(K.row(i)) for i in range(4)]).inv() * K


def dm(Mx):
    return DomainMatrix([[FF.from_sympy(sp.together(x)) for x in Mx.row(i)] for i in range(Mx.rows)],
                        Mx.shape, FF)


Wd, Bd, Id = dm(W), dm(B), dm(sp.eye(4))
tF = FF.from_sympy(t)
P = Wd * (Id - Bd * tF)
R = Id - Wd + (Id + Bd * tF).inv() * (Wd * FF.from_sympy(sp.Integer(2)) - Id)
w_anti = (sp.Rational(1, 10) - sp.Rational(3, 2) * eps) / (2 * (1 + eps))
modes = {'P': [1 - t, w_anti * (1 - t / 2500)],
         'R': [1 / (1 + t), 1 - w_anti + (2 * w_anti - 1) / (1 + t / 2500)]}


def quadratic_factor(Md, roots):
    c = Md.charpoly()
    for r in roots:
        rF = FF.from_sympy(sp.together(r))
        out = [c[0]]
        for ck in c[1:]:
            out.append(ck + out[-1] * rF)
        assert out[-1] == FF.zero
        c = out[:-1]
    return sp.factor(FF.to_sympy(-c[1])), sp.factor(FF.to_sympy(c[2]))


TP, DP = quadratic_factor(P, modes['P'])
TR, DR = quadratic_factor(R, modes['R'])
a0 = eps * (sp.Rational(31, 120) + sp.Rational(217, 200) * eps)
a1 = -sp.Rational(133, 3750000) + sp.Rational(1782779, 10000000) * eps + sp.Rational(969583, 50000000) * eps**2
a2 = sp.Rational(751389, 12500000000) - sp.Rational(15987, 1000000000) * eps - sp.Rational(1582713, 125000000000) * eps**2
Delta = a0 + a1 * t + a2 * t**2
assert sp.cancel((1 - TP + DP) * (1 + eps)**2 - Delta) == 0
theta = sp.factor(sp.together((1 - TR + DR) / Delta))

conditions = {
    'P: 1-D_P': 1 - DP, 'P: 1+T_P+D_P': 1 + TP + DP,
    'R: 1-D_R': 1 - DR, 'R: 1+T_R+D_R': 1 + TR + DR,
    'theta': theta,
    'P antisym: 1-p_a': 1 - modes['P'][1], 'P antisym: 1+p_a': 1 + modes['P'][1],
    'R antisym: 1-r_a': 1 - modes['R'][1], 'R antisym: 1+r_a': 1 + modes['R'][1],
    'a2': a2, '4*a2+a1 (vertex < 2)': 4 * a2 + a1, '-a1': -a1,
    'Delta(eps,2) = a0+2a1+4a2': a0 + 2 * a1 + 4 * a2,
}
side = []
for name, expr in conditions.items():
    num, den = signed_fraction(expr)
    entry = {'name': name, 'numerator': poly2_terms(num), 'denominator': poly2_terms(den)}
    for part, poly in (('numerator', num), ('denominator', den)):
        m, n, arr = bernstein(poly)
        flat = [Fraction(x) for row in arr for x in row]
        entry[part + '_bernstein'] = {'degree_u': m, 'degree_s': n, 'coefficients': arr,
                                      'count': len(flat), 'min': str(min(flat)),
                                      'all_positive': all(x > 0 for x in flat)}
        assert all(x > 0 for x in flat), name + ' ' + part
    side.append(entry)
    print('Bernstein positive on [0,1/10000]x[0,2]:', name)

disc = sp.Poly(sp.expand(a1**2 - 4 * a0 * a2), eps)
iv = [i for i, _ in disc.intervals() if i[1] > 0]
crit = min(iv, key=lambda i: i[0])
crit = disc.refine_root(crit[0], crit[1], eps=sp.Rational(1, 10**15))
assert disc.count_roots(crit[0], crit[1]) == 1 and disc.count_roots(0, EPS1) == 1
thm4 = {
    'description': 'Theorem 4. Box: 0<=epsilon<=1/10000, 0<=t<=2; Bernstein form on the unit square '
                   'after epsilon=u/10000, t=2s. Polynomials are lists [i, j, c] for c*epsilon^i*t^j.',
    'map': {'epsilon': '1/10000 * u', 't': '2 * s'},
    'T_P': str(TP), 'D_P': str(DP), 'T_R': str(TR), 'D_R': str(DR), 'theta': str(theta),
    'Delta': poly2_terms(Delta),
    'side_conditions': side,
    'discriminant_integer_coeffs_ascending': poly1_coeffs(disc),
    'eps_crit_interval': [q(crit[0]), q(crit[1])],
    'eps_crit_decimal': float((crit[0] + crit[1]) / 2),
}
(OUTDIR / 'thm4_classification.json').write_text(json.dumps(thm4, indent=1) + '\n')

fam = sp.Poly(sp.cancel(Delta.subs(eps, t / 100000) / t), t)
tsiv = [i for i, _ in fam.intervals() if i[1] > 0]
ts = min(tsiv, key=lambda i: i[0])
ts = fam.refine_root(ts[0], ts[1], eps=sp.Rational(1, 10**15))
ts = (sp.Rational(5312853253, 10**10), sp.Rational(5312853254, 10**10))
assert fam.count_roots(ts[0], ts[1]) == 1 and fam.count_roots(0, 2) == 1
(OUTDIR / 'thm4_tstar.json').write_text(json.dumps({
    'description': 'Delta(t/100000, t)/t, cleared to integers; t* is its unique root in (0,2].',
    'integer_coeffs_ascending': poly1_coeffs(fam),
    't_star_interval': [q(ts[0]), q(ts[1])], 'check_points': ['0', q(ts[0]), q(ts[1]), '2'],
}, indent=1) + '\n')
print('t* in', [float(x) for x in ts])

# ---------------------------------------------------------------- Theorem 3
n = 12
G1, G2, H1, H2 = [0, 1, 2, 4, 6], [3, 5, 7, 8, 9, 10, 11], {0, 2, 4, 6}, {5, 8, 9, 10}
row = [sp.Rational(x, 144) for x in (33, 11, 1, 11, 21, 11, 1, 11, 21, 11, 1, 11)]
Ac = sp.Matrix(n, n, lambda i, j: row[(j - i) % n])
Bc = Ac * Ac


def k0(i, j):
    if i == j:
        return sp.Integer(1)
    for G, Hh in ((G1, H1), (G2, H2)):
        if i in G and j in G:
            return sp.Rational(99, 100) if (i in Hh and j in Hh) else sp.Rational(1, 100)
    return sp.Integer(0)


Kc = (1 - t**2) * sp.Matrix(n, n, k0) + t**2 * sp.ones(n, n)
Dc = sp.diag(*[sp.expand(sum(Kc.row(i))) for i in range(n)])
M = DomainMatrix.from_Matrix(Dc - Kc + t * Kc * Bc).convert_to(sp.QQ[t])
Nc = sp.Poly(sp.QQ[t].to_sympy(M.det()), t)
assert Nc.degree() == 36 and Nc.eval(0) == 0
pc = sp.Poly(sp.cancel(Nc.as_expr() / t**2), t)
L, U = sp.Rational(19157074941, 10**14), sp.Rational(19157074942, 10**14)
assert pc.count_roots(L, U) == 1 and pc.count_roots(0, 2) == 1 and pc.eval(0) < 0
(OUTDIR / 'thm3_circulant.json').write_text(json.dumps({
    'description': 'p(t) = N_c(t)/t^2 with N_c(t) = det(D - K + t K B_c) along epsilon = t^2, '
                   'cleared to integer coefficients (a positive multiple). t_c is its unique root in (0,2].',
    'integer_coeffs_ascending': poly1_coeffs(pc),
    't_c_interval': [q(L), q(U)], 'check_points': ['0', q(L), q(U), '1', '2'],
}, indent=1) + '\n')
print('t_c in', [float(L), float(U)], 'degree of p:', pc.degree())
print('certificate files written to', OUTDIR.name)
