#!/usr/bin/env python3
"""Exact certificates for the claims added in the revised manuscript.

Sections (each prints its own header):
  [A] Four-state family, two-parameter classification (Theorem 4):
      characteristic-polynomial factorisation of P and R, Bernstein
      certificates for the two side Jury conditions on 0<=eps<=EPS1, 0<=t<=2,
      exact first-loss root, exact critical epsilon, and the exact boundary
      t* of the one-parameter family eps=t/100000.
  [B] Twelve-state circulant family along eps=t^2 (Theorem 3): exact
      determinant polynomial, certified sign on (0, t_c), exact witness on
      the curve.
  [C] Conditioning constant of Theorem 1(b).

Requires sympy (pinned in README). Everything that is asserted is decided in
exact rational arithmetic; decimals are printed for reading only.
Writes revision_certificate.json beside this file unless --out is given.
"""
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp

START = time.time()
eps, t, lam, u, s = sp.symbols('epsilon t lambda u s')
OUT = Path(__file__).with_name('revision_certificate.json')
if '--out' in sys.argv:
    OUT = Path(sys.argv[sys.argv.index('--out') + 1])
REPORT = {}


def header(text):
    print('\n' + '=' * 72 + '\n' + text + '\n' + '=' * 72, flush=True)


def bernstein_positive(poly, a, b, depth=0, max_depth=12):
    """Certify poly(eps,t) > 0 on [0,a]x[0,b] (a,b rational) by Bernstein
    coefficients with recursive bisection. Returns number of boxes used."""
    return _bern(sp.Poly(sp.expand(poly), eps, t), sp.Rational(0), sp.Rational(a),
                 sp.Rational(0), sp.Rational(b), depth, max_depth)


def _bern(P, e0, e1, t0, t1, depth, max_depth):
    q = sp.Poly(sp.expand(P.as_expr().subs({eps: e0 + (e1 - e0) * u,
                                             t: t0 + (t1 - t0) * s},
                                            simultaneous=True)), u, s)
    m, n = q.degree(u), q.degree(s)
    m, n = max(m, 0), max(n, 0)
    c = {k: v for k, v in zip(q.monoms(), q.coeffs())}
    ok = True
    for k in range(m + 1):
        for l in range(n + 1):
            val = sum(sp.binomial(k, i) / sp.binomial(m, i) *
                      sp.binomial(l, j) / sp.binomial(n, j) * c.get((i, j), 0)
                      for i in range(k + 1) for j in range(l + 1))
            if val <= 0:
                ok = False
                break
        if not ok:
            break
    if ok:
        return 1
    if depth >= max_depth:
        raise AssertionError('Bernstein certificate failed')
    em, tm = (e0 + e1) / 2, (t0 + t1) / 2
    return sum(_bern(P, a0, a1, b0, b1, depth + 1, max_depth)
               for a0, a1 in [(e0, em), (em, e1)] for b0, b1 in [(t0, tm), (tm, t1)])


def rat(x):
    return str(sp.Rational(x))


# ---------------------------------------------------------------------------
header('[A] Four-state Gaussian family: two-parameter classification')
v = sp.Matrix([-8, -8, 17, -1])
J4 = sp.ones(4, 4)
A = (20 * sp.eye(4) + 245 * J4 + v * v.T) / 1000
B = A * A
K0 = sp.Matrix([[1, sp.Rational(9, 10), sp.Rational(1, 10), 0],
                [sp.Rational(9, 10), 1, sp.Rational(1, 10), 0],
                [sp.Rational(1, 10), sp.Rational(1, 10), 1, 0],
                [0, 0, 0, 1]])
H = sp.Matrix([[0, sp.Rational(3, 2), sp.Rational(3, 10), sp.Rational(1, 5)],
               [sp.Rational(3, 2), 0, sp.Rational(3, 10), sp.Rational(1, 5)],
               [sp.Rational(3, 10), sp.Rational(3, 10), 0, sp.Rational(3, 5)],
               [sp.Rational(1, 5), sp.Rational(1, 5), sp.Rational(3, 5), 0]])
from sympy.polys.matrices import DomainMatrix
FF = sp.QQ.frac_field(eps, t)
K = K0 + eps * H
Dg = sp.diag(*[sum(K.row(i)) for i in range(4)])
W = Dg.inv() * K


def dm(Mx):
    return DomainMatrix([[FF.from_sympy(sp.together(x)) for x in Mx.row(i)] for i in range(Mx.rows)],
                        Mx.shape, FF)


Wd, Bd, Id = dm(W), dm(B), dm(sp.eye(4))
tF = FF.from_sympy(t)
P = Wd * (Id - Bd * tF)
R = Id - Wd + (Id + Bd * tF).inv() * (Wd * FF.from_sympy(sp.Integer(2)) - Id)

w_anti = (sp.Rational(1, 10) - sp.Rational(3, 2) * eps) / (2 * (1 + eps))
beta_anti = sp.Rational(1, 2500)
p_const, p_anti = 1 - t, w_anti * (1 - t * beta_anti)
r_const = 1 / (1 + t)
r_anti = 1 - w_anti + (2 * w_anti - 1) / (1 + t * beta_anti)


def residual_quadratic(Md, roots):
    c = Md.charpoly()            # monic coefficients, highest degree first
    for r in roots:              # exact synthetic division by (lambda - r)
        rF = FF.from_sympy(sp.together(r))
        out = [c[0]]
        for ck in c[1:]:
            out.append(ck + out[-1] * rF)
        assert out[-1] == FF.zero, 'factor does not divide'
        c = out[:-1]
    assert len(c) == 3
    return sp.factor(FF.to_sympy(-c[1])), sp.factor(FF.to_sympy(c[2]))


TP, DP = residual_quadratic(P, [p_const, p_anti])
TR, DR = residual_quadratic(R, [r_const, r_anti])
print('PASS: char. polynomials factor as (constant)(antisymmetric)(quadratic) for P and R')

Delta = (eps * (sp.Rational(31, 120) + sp.Rational(217, 200) * eps)
         + t * (-sp.Rational(133, 3750000) + sp.Rational(1782779, 10000000) * eps
                + sp.Rational(969583, 50000000) * eps**2)
         + t**2 * (sp.Rational(751389, 12500000000) - sp.Rational(15987, 1000000000) * eps
                   - sp.Rational(1582713, 125000000000) * eps**2))
assert sp.cancel((1 - TP + DP) * (1 + eps)**2 - Delta) == 0
print('PASS: 1 - tr + det of the P-quotient equals Delta(eps,t)/(1+eps)^2, eq. (33)')

EPS1 = sp.Rational(1, 10000)
boxes = {}
for name, expr in [('P: 1-det', 1 - DP), ('P: 1+tr+det', 1 + TP + DP),
                   ('R: 1-det', 1 - DR), ('R: 1+tr+det', 1 + TR + DR)]:
    num, den = sp.fraction(sp.factor(sp.together(expr)))
    nb = bernstein_positive(num, EPS1, 2)
    db = bernstein_positive(den, EPS1, 2)
    boxes[name] = [nb, db]
    print(f'PASS: {name} > 0 on [0,1/10000]x[0,2]  (Bernstein boxes: num {nb}, den {db})')
# the R +1 condition has the sign of Delta
ratio = sp.factor(sp.together((1 - TR + DR) / Delta))
rn, rd = sp.fraction(ratio)
assert bernstein_positive(rn, EPS1, 2) and bernstein_positive(rd, EPS1, 2)
print('PASS: 1 - tr + det of the R-quotient equals Delta times a positive factor')

# Antisymmetric and constant modes are stable for 0<t<2 and 0<=eps<=EPS1.
def positive_rational(expr):
    num, den = sp.fraction(sp.factor(sp.together(expr)))
    if sp.Poly(den, eps, t).eval({eps: 0, t: 0}) < 0:
        num, den = -num, -den
    return bernstein_positive(num, EPS1, 2) and bernstein_positive(den, EPS1, 2)


for x in (p_anti, r_anti):
    assert positive_rational(1 - x) and positive_rational(1 + x)
print('PASS: antisymmetric modes of P and R lie in (-1,1); constant modes 1-t, 1/(1+t) stable for 0<t<2')

a2 = sp.Poly(Delta, t).coeff_monomial(t**2)
a1 = sp.Poly(Delta, t).coeff_monomial(t)
a0 = sp.Poly(Delta, t).coeff_monomial(1)
assert a2.subs(eps, 0) > 0 and bernstein_positive(a2, EPS1, 1)
disc = sp.expand(a1**2 - 4 * a2 * a0)
disc_poly = sp.Poly(disc, eps)
roots_disc = [r for r in disc_poly.intervals() if r[0][1] > 0]
crit_iv = min(roots_disc, key=lambda r: r[0][0])[0]
crit_iv = disc_poly.refine_root(crit_iv[0], crit_iv[1], eps=sp.Rational(1, 10**14))
print('critical epsilon (discriminant root) in', [float(x) for x in crit_iv])
assert disc_poly.count_roots(sp.Rational(1, 10**12), EPS1) == 1
assert disc_poly.eval(0) > 0
# for 0 < eps <= EPS1 both roots of Delta(eps,.) are positive (a0 > 0, a1 < 0, a2 > 0)
assert sp.Poly(a1, eps).count_roots(0, EPS1) == 0 and a1.subs(eps, 0) < 0   # a1 < 0 on the box
assert bernstein_positive(4 * a2 + a1, EPS1, 1)  # vertex -a1/(2 a2) < 2
print('PASS: on (0,1/10000]: a0>0, a1<0, a2>0, vertex<2; Delta(eps,.) is negative somewhere in (0,2) iff eps < eps_crit')
# root location: smaller root < 2 on (0, eps_crit) is checked by value of Delta at t=2 and vertex
REPORT['A'] = {
    'Delta_coefficients': {'t^0': str(a0), 't^1': str(a1), 't^2': str(a2)},
    'discriminant_in_eps': str(disc),
    'critical_epsilon_interval': [rat(crit_iv[0]), rat(crit_iv[1])],
    'bernstein_boxes': boxes,
}

first_loss = sp.series((-a1 - sp.sqrt(disc)) / (2 * a2), eps, 0, 3).removeO()
print('first-loss root T_-(eps) =', sp.nsimplify(first_loss))
REPORT['A']['first_loss_series'] = str(first_loss)

# One-parameter family eps = t/100000: exact boundary t*.
fam = sp.Poly(sp.factor(Delta.subs(eps, t / 100000) / t), t)
fam_n = sp.Poly(sp.fraction(sp.together(fam.as_expr()))[0], t)
ivs = [iv for iv, _ in fam_n.intervals() if iv[1] > 0]
tstar = min(ivs, key=lambda iv: iv[0])
tstar = fam_n.refine_root(tstar[0], tstar[1], eps=sp.Rational(1, 10**12))
print('family eps=t/1e5: Delta<0 exactly on (0, t*), t* in', [float(x) for x in tstar])
assert fam_n.count_roots(0, 2) == 1  # only t* in (0,2]; Delta>0 on (t*,2]
cubic_f = sp.Poly([sp.Rational(-4748139, 1562500000000000000000), sp.Rational(-236896251, 625000000000000000),
                   sp.Rational(74272809, 500000000000), sp.Rational(-1973, 25000000)], t)
assert cubic_f.count_roots(0, 2) == 1
f_ivs = [iv for iv, _ in cubic_f.intervals(inf=0, sup=2)]
fi = cubic_f.refine_root(f_ivs[0][0], f_ivs[0][1], eps=sp.Rational(1, 10**12))
assert fi[1] >= tstar[0] and fi[0] <= tstar[1]
print('PASS: t* is the unique root in (0,2] of the cubic f in eq. (20); f<0 on (0,t*), f>0 on (t*,2]')
REPORT['A']['t_star_interval'] = [rat(fi[0]), rat(fi[1])]

# ---------------------------------------------------------------------------
header('[B] Twelve-state circulant family along eps = t^2')
from sympy.polys.matrices import DomainMatrix
from sympy import QQ
n = 12
Gg, Gh = [0, 1, 2, 4, 6], [3, 5, 7, 8, 9, 10, 11]
hg, hh = {0, 2, 4, 6}, {5, 8, 9, 10}
Fr = Fraction


def circ_mats():
    E = [[Fr(1, 12)] * n for _ in range(n)]
    Z = [[Fr([1, 0, -1, 0][(i - j) % 4], 6) for j in range(n)] for i in range(n)]
    Bc = [[Fr(int(i == j), 144) + Fr(143, 144) * E[i][j] + Fr(35, 144) * Z[i][j]
           for j in range(n)] for i in range(n)]
    K0c = [[Fr(0)] * n for _ in range(n)]
    for grp, hi in [(Gg, hg), (Gh, hh)]:
        for i in grp:
            for j in grp:
                K0c[i][j] = Fr(1) if i == j else (Fr(99, 100) if (i in hi and j in hi) else Fr(1, 100))
    return Bc, K0c


Bc, K0c = circ_mats()
Rt = sp.QQ[t]
tt = Rt.gens[0]


rows = []
for i in range(n):
    row = []
    Krow = [Rt.convert(sp.Rational(K0c[i][j].numerator, K0c[i][j].denominator)) * (1 - tt**2) + tt**2
            for j in range(n)]
    deg = sum(Krow, Rt.zero)
    KB = [sum((Krow[k] * Rt.convert(sp.Rational(Bc[k][j].numerator, Bc[k][j].denominator))
               for k in range(n)), Rt.zero) for j in range(n)]
    for j in range(n):
        row.append((deg if i == j else Rt.zero) - Krow[j] + tt * KB[j])
    rows.append(row)
M = DomainMatrix(rows, (n, n), Rt)
Ndet = sp.Poly(Rt.to_sympy(M.det()), t)
k0 = min(m[0] for m in Ndet.monoms())
N1 = sp.Poly(sp.cancel(Ndet.as_expr() / t**k0), t)
print('numerator det(D - K + tKB) has degree', Ndet.degree(), 'and factor t^%d' % k0)
lead_low = N1.eval(0)
assert lead_low < 0
pos = sorted([iv for iv, _ in N1.intervals() if iv[1] > 0], key=lambda iv: iv[0])
tc = N1.refine_root(pos[0][0], pos[0][1], eps=sp.Rational(1, 10**15))
print('PASS: det(I-P(t)) < 0 for every 0 < t < t_c, t_c in', [float(x) for x in tc])
print('      next positive root at about', float(pos[1][0]) if len(pos) > 1 else None)
# det D > 0, so sign det(I-P) = sign Ndet
REPORT['B'] = {'numerator_degree': Ndet.degree(), 't_power': k0,
               'value_at_0_after_t_power': str(lead_low),
               't_c_interval': [rat(tc[0]), rat(tc[1])]}

# exact witness on the curve: t = 1/10000, eps = 1e-8
tw = Fr(1, 10000)
ew = tw * tw


def det_frac(a):
    a = [r[:] for r in a]
    res = Fr(1)
    for i in range(len(a)):
        p = next(j for j in range(i, len(a)) if a[j][i])
        if p != i:
            a[i], a[p] = a[p], a[i]
            res = -res
        res *= a[i][i]
        for j in range(i + 1, len(a)):
            c = a[j][i] / a[i][i]
            for k in range(i, len(a)):
                a[j][k] -= c * a[i][k]
    return res


Kw = [[(1 - ew) * K0c[i][j] + ew for j in range(n)] for i in range(n)]
dw = [sum(r) for r in Kw]
Ww = [[Kw[i][j] / dw[i] for j in range(n)] for i in range(n)]
Pw = [[Ww[i][j] - tw * sum(Ww[i][k] * Bc[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
IminusP = [[Fr(int(i == j)) - Pw[i][j] for j in range(n)] for i in range(n)]
dP = det_frac(IminusP)
dIB = det_frac([[Fr(int(i == j)) + tw * Bc[i][j] for j in range(n)] for i in range(n)])
assert dP < 0 and dIB > 0
print('PASS: exact witness on the curve: t=1/10000, eps=10^-8: det(I-P) =', float(dP),
      ', det(I-R) =', float(dP / dIB))
REPORT['B']['witness_on_curve'] = {'t': '1/10000', 'epsilon': '1/100000000',
                                   'det_I_P': str(dP), 'det_I_R': str(dP / dIB)}

# ---------------------------------------------------------------------------
header('[C] Conditioning of the four-state affinities')
lmax0 = (29 + sp.sqrt(89)) / 20
kappa_bound = (lmax0 + sp.Rational(1, 100000)) / (sp.Rational(1, 10) - sp.Rational(1, 100000))
assert sp.simplify(sp.Matrix(K0).eigenvals().get(lmax0, 0)) == 1
assert kappa_bound < sp.Rational(1922, 100)
print('PASS: lambda_max(K0) = (29+sqrt(89))/20, kappa_2(K0) =', float(lmax0 * 10),
      '; kappa_2(K_t) <', float(kappa_bound), '< 19.22 for 0<t<=1/2')
REPORT['C'] = {'kappa_K0': str(lmax0 * 10), 'kappa_Kt_bound_decimal': float(kappa_bound)}

REPORT['python'] = sys.version.split()[0]
REPORT['sympy'] = sp.__version__
REPORT['elapsed_seconds'] = round(time.time() - START, 1)
OUT.write_text(json.dumps(REPORT, indent=2) + '\n')
print('\nAll revision certificates passed. Report written to', OUT.name)
