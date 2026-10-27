#!/usr/bin/env python3
"""Certificate-only checker for Theorems 3 and 4 (Python standard library only).

Reads certificates/*.json (written by emit_certificates.py) and, without any
symbolic derivation, re-verifies:
  * every Bernstein array is the Bernstein form of its stated polynomial under
    epsilon = u/10000, t = 2s, and every coefficient is strictly positive;
  * the discriminant and the t* polynomial follow from the stated Delta;
  * Sturm sign-variation counts: one root in each isolating interval, and the
    stated root is the only one in (0, 1/10000] (discriminant) or (0, 2].
Prints SHA-256 hashes of the files it read. Exit status 0 means every check passed.
"""
import hashlib
import json
import sys
from fractions import Fraction as F
from math import comb
from pathlib import Path

D = Path(__file__).resolve().parent / 'certificates'
if '--dir' in sys.argv:
    D = Path(sys.argv[sys.argv.index('--dir') + 1])
EPS1, TMAX = F(1, 10000), F(2)


def load(name):
    raw = (D / name).read_bytes()
    print(f'{hashlib.sha256(raw).hexdigest()}  {name}')
    return json.loads(raw)


def padd(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v}


def pmul(a, b):
    out = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            out[(i + k, j + l)] = out.get((i + k, j + l), 0) + x * y
    return {k: v for k, v in out.items() if v}


def terms(lst):
    return {(i, j): F(c) for i, j, c in lst}


def bernstein(p):
    m = max([i for i, _ in p] + [0])
    n = max([j for _, j in p] + [0])
    c = {(i, j): x * EPS1**i * TMAX**j for (i, j), x in p.items()}
    return [[sum(F(comb(k, i), comb(m, i)) * F(comb(l, j), comb(n, j)) * c.get((i, j), 0)
                 for i in range(k + 1) for j in range(l + 1)) for l in range(n + 1)]
            for k in range(m + 1)]


def strip(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def prem(a, b):
    a = list(a)
    while len(a) >= len(b) and a:
        c = a[-1] / b[-1]
        shift = len(a) - len(b)
        for i, x in enumerate(b):
            a[i + shift] -= c * x
        a = strip(a)
    return a


def sturm(p):
    seq = [strip(p), strip([i * c for i, c in enumerate(p)][1:])]
    while len(seq[-1]) > 1:
        r = prem(seq[-2], seq[-1])
        if not r:
            break  # seq[-1] is gcd(p, p'); counts below are of distinct real roots
        seq.append([-x for x in r])
    return seq


def variations(seq, x):
    vals = [sum(c * x**i for i, c in enumerate(p)) for p in seq]
    vals = [v for v in vals if v != 0]
    return sum(1 for a, b in zip(vals, vals[1:]) if (a > 0) != (b > 0))


def ev(p, x):
    return sum(c * x**i for i, c in enumerate(p))


def roots_between(seq, a, b):
    return variations(seq, a) - variations(seq, b)


ok = True


def check(cond, msg):
    global ok
    print(('PASS: ' if cond else 'FAIL: ') + msg)
    ok = ok and cond


# ------------------------------------------------------------ Theorem 4
c4 = load('thm4_classification.json')
for e in c4['side_conditions']:
    for part in ('numerator', 'denominator'):
        arr = bernstein(terms(e[part]))
        stored = [[F(x) for x in row] for row in e[part + '_bernstein']['coefficients']]
        check(arr == stored and all(x > 0 for row in arr for x in row),
              f"{e['name']} {part}: {sum(len(r) for r in arr)} Bernstein coefficients, all > 0, "
              f"min {float(min(x for r in arr for x in r)):.3e}")
Delta = terms(c4['Delta'])
a = [{(i, 0): x for (i, j), x in Delta.items() if j == k} for k in range(3)]
disc = padd(pmul(a[1], a[1]), {k: -4 * v for k, v in pmul(a[0], a[2]).items()})
dpoly = [F(int(x)) for x in c4['discriminant_integer_coeffs_ascending']]
dfrom = [disc.get((i, 0), F(0)) for i in range(len(dpoly))]
scale = dpoly[0] / dfrom[0]
check(scale > 0 and all(x * scale == y for x, y in zip(dfrom, dpoly)),
      'discriminant = a1^2 - 4 a0 a2 (positive multiple)')
S = sturm(dpoly)
L, U = (F(x) for x in c4['eps_crit_interval'])
check(len(S[-1]) == 1 and ev(dpoly, L) * ev(dpoly, U) < 0, 'discriminant square-free; sign change across interval')
check(roots_between(S, L, U) == 1 and roots_between(S, F(0), EPS1) == 1 and ev(dpoly, F(0)) > 0,
      f'eps_crit: Sturm counts V(0)={variations(S, F(0))}, V(L)={variations(S, L)}, '
      f'V(U)={variations(S, U)}, V(1e-4)={variations(S, EPS1)}; unique root in ({L}, {U})')

ts = load('thm4_tstar.json')
fp = [F(int(x)) for x in ts['integer_coeffs_ascending']]
fam = {}
for (i, j), x in Delta.items():
    fam[i + j - 1] = fam.get(i + j - 1, 0) + x / F(100000)**i
fl = [fam.get(k, F(0)) for k in range(len(fp))]
sc = fp[0] / fl[0]
check(sc > 0 and all(x * sc == y for x, y in zip(fl, fp)), 't* polynomial = Delta(t/100000, t)/t (positive multiple)')
S = sturm(fp)
L, U = (F(x) for x in ts['t_star_interval'])
check(len(S[-1]) == 1 and ev(fp, L) < 0 < ev(fp, U), 't* polynomial square-free; sign change across interval')
check(roots_between(S, L, U) == 1 and roots_between(S, F(0), F(2)) == 1 and ev(fp, F(0)) < 0,
      f'Delta(t/1e5,t) < 0 on (0,t*), > 0 on (t*,2]: V(0)={variations(S, F(0))}, V(L)={variations(S, L)}, '
      f'V(U)={variations(S, U)}, V(2)={variations(S, F(2))}; t* in ({float(L)}, {float(U)})')

# ------------------------------------------------------------ Theorem 3
c3 = load('thm3_circulant.json')
pp = [F(int(x)) for x in c3['integer_coeffs_ascending']]
S = sturm(pp)
L, U = (F(x) for x in c3['t_c_interval'])
print(f'      gcd(p, p\') has degree {len(S[-1]) - 1} (repeated factors; counts are of distinct roots)')
check(ev(pp, L) < 0 < ev(pp, U), 'p changes sign across the isolating interval: p(L) < 0 < p(U)')
check(len(pp) == 35 and ev(pp, F(0)) < 0 and roots_between(S, L, U) == 1 and roots_between(S, F(0), F(2)) == 1,
      f'p(t)=N_c/t^2, degree {len(pp) - 1}: p(0)<0; V(0)={variations(S, F(0))}, V(L)={variations(S, L)}, '
      f'V(U)={variations(S, U)}, V(1)={variations(S, F(1))}, V(2)={variations(S, F(2))}; '
      f't_c in ({float(L)}, {float(U)}), the only root in (0,2]')

print('ALL CERTIFICATE CHECKS PASSED' if ok else 'CERTIFICATE CHECK FAILED')
sys.exit(0 if ok else 1)
