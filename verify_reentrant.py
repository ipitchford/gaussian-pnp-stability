#!/usr/bin/env python3
"""Exact algebra checks for Theorem 5 (Section 10) of gaussian_pnp_stability.md:
the three-state reentrant family with unequal diagonal entries.

Standard-library only.  Polynomial identities are compared coefficient by
coefficient over Fraction; no floating-point eigenvalue is used as a certificate.
Run from any directory: python3 verify_reentrant.py [--out FILE]
(default output: reentrant_verification.json beside this script).
"""
from fractions import Fraction as Q
from itertools import permutations
from pathlib import Path
import json
import platform
import time


class Poly:
    """Sparse Q[epsilon,t,lambda], deliberately small and transparent."""
    def __init__(self, value=0):
        if isinstance(value, Poly):
            self.c = value.c.copy()
        elif isinstance(value, dict):
            self.c = {k: Q(v) for k, v in value.items() if v}
        else:
            self.c = {} if not value else {(0, 0, 0): Q(value)}

    def __add__(self, other):
        other = Poly(other)
        out = self.c.copy()
        for k, v in other.c.items():
            out[k] = out.get(k, Q(0)) + v
        return Poly(out)

    __radd__ = __add__

    def __neg__(self):
        return Poly({k: -v for k, v in self.c.items()})

    def __sub__(self, other):
        return self + (-Poly(other))

    def __rsub__(self, other):
        return Poly(other) - self

    def __mul__(self, other):
        other = Poly(other)
        out = {}
        for a, u in self.c.items():
            for b, v in other.c.items():
                k = tuple(x + y for x, y in zip(a, b))
                out[k] = out.get(k, Q(0)) + u * v
        return Poly(out)

    __rmul__ = __mul__

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        out = Poly(1)
        for _ in range(n):
            out *= self
        return out

    def at(self, epsilon=0, t=0, lam=0):
        x = [Q(epsilon), Q(t), Q(lam)]
        return sum((v * x[0] ** k[0] * x[1] ** k[1] * x[2] ** k[2]
                    for k, v in self.c.items()), Q(0))

    def serialise(self):
        return [{"powers_epsilon_t_lambda": list(k), "coefficient": str(v)}
                for k, v in sorted(self.c.items())]


def add(a, b):
    return [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(s, a):
    return [[s * x for x in row] for row in a]


def mul(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(len(b))), Poly(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def det(a):
    n = len(a)
    ans = Poly(0)
    for p in permutations(range(n)):
        inversions = sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        term = Poly((-1) ** inversions)
        for i in range(n):
            term *= a[i][p[i]]
        ans += term
    return ans


def eye(n):
    return [[Poly(i == j) for j in range(n)] for i in range(n)]


def const(a):
    return [[Poly(x) for x in row] for row in a]


def exact_equal(a, b, name, checks):
    residual = a - b
    assert not residual.c, (name, residual.serialise())
    checks.append({"name": name, "status": "passed", "residual_coefficients": 0})


def main():
    started = time.monotonic()
    eps = Poly({(1, 0, 0): 1})
    t = Poly({(0, 1, 0): 1})
    lam = Poly({(0, 0, 1): 1})
    tau = 1 - eps
    I = eye(3)
    B = const([[Q(x, 180) for x in row]
               for row in [[68, 71, 41], [71, 86, 23], [41, 23, 116]]])
    S = const([[1, 0, 0], [0, Q(19, 20), Q(1, 20)],
               [0, Q(19, 20), Q(1, 20)]])
    pi = [Q(1, 2), Q(19, 40), Q(1, 40)]
    Epi = const([pi, pi, pi])
    W = add(scale(tau * Q(1, 2), add(I, S)), scale(eps, Epi))
    P = mul(W, add(I, scale(-t, B)))
    checks = []
    for i in range(3):
        exact_equal(sum(W[i], Poly(0)), Poly(1), f"W row {i+1} stochastic", checks)
        exact_equal(sum(B[i], Poly(0)), Poly(1), f"B row {i+1} stochastic", checks)
        exact_equal(sum((pi[j] * W[j][i] for j in range(3)), Poly(0)),
                    Poly(pi[i]), f"pi stationary coordinate {i+1}", checks)
        for j in range(3):
            exact_equal(pi[i] * W[i][j], pi[j] * W[j][i],
                        f"detailed balance ({i+1},{j+1})", checks)

    exact_equal(det(add(scale(lam, I), scale(-1, B))),
                (lam-1)*(lam-Q(29,60))*(lam-Q(1,60)), "B characteristic polynomial", checks)
    exact_equal(det(add(scale(lam, I), scale(-1, W))),
                (lam-1)*(lam-tau)*(lam-tau*Q(1,2)),
                "W characteristic polynomial", checks)

    h, d, c = Q(1, 2), Q(29, 3600), -Q(1, 120)
    sq, dq, u = tau * Q(3, 2), tau**2 * Q(1, 2), tau * Q(59, 240)
    traceP = sq-t*u
    determinantP = dq*(1-h*t+d*t*t)
    F = eps*(1+eps)*Q(1,2) + tau*(60*eps-1)*t*Q(1,240) + tau*tau*t*t*Q(29,7200)
    exact_equal(det(add(scale(lam, I), scale(-1, P))),
                (lam-1+t)*(lam*lam-traceP*lam+determinantP),
                "full P characteristic polynomial", checks)
    exact_equal(det(add(I, scale(-1, P))), t*F, "det(I-P)=t F", checks)
    exact_equal(1-traceP+determinantP, F, "P quotient Jury value at +1", checks)

    G = 1+h*t+d*t*t
    traceN = sq+t*(h*(1+sq)-2*u)+d*t*t*(2-sq)
    determinantN = dq+t*(h*(sq-dq)-u)+d*t*t*(1-sq+dq)
    # R = (I+tB)^(-1) [tB+(I-tB)W].  Avoid rational functions:
    # det(lambda (I+tB)-[...]) is a polynomial identity.
    N = add(scale(t, B), mul(add(I, scale(-t, B)), W))
    pencil = add(scale(lam, add(I, scale(t, B))), scale(-1, N))
    exact_equal(det(pencil), ((1+t)*lam-1)*(G*lam*lam-traceN*lam+determinantN),
                "full R generalised characteristic polynomial", checks)
    exact_equal(G-traceN+determinantN, F, "R quotient Jury value at +1", checks)
    exact_equal(det(add(I, scale(t, B))), (1+t)*G, "positive R denominator", checks)

    # The quadratic discriminant is exact and its root formula is algebraic.
    a = tau*tau*Q(29,7200)
    b = tau*(60*eps-1)*Q(1,240)
    cc = eps*(1+eps)*Q(1,2)
    discr = 3136*eps*eps-584*eps+1
    exact_equal(b*b-4*a*cc, tau*tau*discr*Q(1,57600),
                "phase-boundary discriminant", checks)

    concrete_eps, concrete_t = Q(1,1000), Q(1,4)
    fval = F.at(concrete_eps, concrete_t)
    pdet = concrete_t*fval
    rdenom = ((1+t)*G).at(concrete_eps, concrete_t)
    rdet = pdet/rdenom
    assert pdet < 0 and rdet < 0
    numeric_matrix = lambda M: [[str(x.at(concrete_eps, concrete_t)) for x in row] for row in M]

    # B is written in orthogonal projections.  This verifies A^2=B
    # without evaluating square roots: E and V are orthogonal projectors.
    E = const([[Q(1,3)]*3 for _ in range(3)])
    v = [1,2,-3]
    V = const([[Q(v[i]*v[j],14) for j in range(3)] for i in range(3)])
    T = add(add(I, scale(-1,E)), scale(-1,V))
    for name, M in [("consensus projector", E), ("v projector", V), ("third projector", T)]:
        square = mul(M,M)
        for i in range(3):
            for j in range(3):
                exact_equal(square[i][j], M[i][j], f"{name} idempotent ({i},{j})", checks)
    decomposition = add(add(E, scale(Q(29,60),V)), scale(Q(1,60),T))
    for i in range(3):
        for j in range(3):
            exact_equal(decomposition[i][j], B[i][j], f"B spectral decomposition ({i},{j})", checks)

    report = {
        "status":"All listed exact checks passed; unrestricted theorems have prose proofs, not proof-assistant verification.",
        "python":platform.python_version(), "elapsed_seconds":time.monotonic()-started,
        "checks":checks,
        "concrete":{"epsilon":str(concrete_eps),"t":str(concrete_t),"B":numeric_matrix(B),
                    "W":numeric_matrix(W),"pi":[str(x) for x in pi],
                    "det_I_minus_P":str(pdet),"det_I_minus_R":str(rdet)},
        "phase_polynomial_F":F.serialise(),
        "discriminant_polynomial":discr.serialise(),
        "scope_limits":["No formal proof assistant was used.",
                        "Universal sign inequalities are proved in Section 10 of gaussian_pnp_stability.md; no parameter grid is used.",
                        "Novelty is assessed separately against the recorded literature corpus."]
    }
    import sys
    destination = (Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv
                   else Path(__file__).with_name("reentrant_verification.json"))
    destination.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status":"passed","checks":len(checks),
                      "det_I_minus_P":str(pdet),"det_I_minus_R":str(rdet),
                      "report":destination.name,"elapsed_seconds":report["elapsed_seconds"]},indent=2))


if __name__ == "__main__":
    main()
