#!/usr/bin/env python3
"""Growth rates, affinity sizes and a symmetric-denoiser sanity check (Section 12).
High-precision floating point (mpmath, 50 digits). Illustrative, not a proof."""
import math
import mpmath as mp
mp.mp.dps = 50
v = mp.matrix([-8, -8, 17, -1])
I4 = mp.eye(4)
A = (20 * I4 + 245 * mp.ones(4, 4) + v * v.T) / 1000
B = A * A
K0 = mp.matrix([[1, .9, .1, 0], [.9, 1, .1, 0], [.1, .1, 1, 0], [0, 0, 0, 1]])
H = mp.matrix([[0, 1.5, .3, .2], [1.5, 0, .3, .2], [.3, .3, 0, .6], [.2, .2, .6, 0]])
K0 = mp.matrix([[mp.mpf(str(x)) for x in row] for row in K0.tolist()])
H = mp.matrix([[mp.mpf(str(x)) for x in row] for row in H.tolist()])


def rho(M):
    return max(abs(z) for z in mp.eig(M)[0])


print('=' * 60 + '\n[1] Four-state family, eps = t/100000\n' + '=' * 60)
print(' t        rho(P)-1     rho(R)-1     iters for 10x   smallest affinity')
for t in ['0.000001', '0.01', '0.1', '0.25', '0.5']:
    t = mp.mpf(t); e = t / 100000
    K = K0 + e * H
    W = mp.matrix(4, 4)
    for i in range(4):
        d = sum(K[i, j] for j in range(4))
        for j in range(4):
            W[i, j] = K[i, j] / d
    P = W * (I4 - t * B)
    R = I4 - W + mp.inverse(I4 + t * B) * (2 * W - I4)
    gp, gr = rho(P) - 1, rho(R) - 1
    print(f' {float(t):<8g} {float(gp):.3e}   {float(gr):.3e}   {math.log(10)/math.log1p(float(gp)):.3e}       {float(e/5):.1e}')
    if float(t) == 0.25:
        Dh = mp.diag([1 / mp.sqrt(sum(K[i, j] for j in range(4))) for i in range(4)])
        S = Dh * K * Dh
        print('   symmetric normalisation at t=1/4: rho(S(I-tB)) - 1 =', mp.nstr(rho(S * (I4 - t * B)) - 1, 5))

print('=' * 60 + '\n[2] Twelve-state circulant family, eps = t^2\n' + '=' * 60)
n = 12
Gg, Gh, hg, hh = [0, 1, 2, 4, 6], [3, 5, 7, 8, 9, 10, 11], {0, 2, 4, 6}, {5, 8, 9, 10}
Bc = mp.matrix(n, n)
for i in range(n):
    for j in range(n):
        Bc[i, j] = mp.mpf(int(i == j)) / 144 + mp.mpf(143) / 144 / 12 + mp.mpf(35) / 144 * [1, 0, -1, 0][(i - j) % 4] / 6
K0c = mp.matrix(n, n)
for grp, hi in [(Gg, hg), (Gh, hh)]:
    for i in grp:
        for j in grp:
            K0c[i, j] = 1 if i == j else (mp.mpf(99) / 100 if (i in hi and j in hi) else mp.mpf(1) / 100)
for t in ['0.00005', '0.0001', '0.00015', '0.0005', '0.01']:
    t = mp.mpf(t); e = t * t
    K = (1 - e) * K0c + e * mp.ones(n, n)
    W = mp.matrix(n, n)
    for i in range(n):
        d = sum(K[i, j] for j in range(n))
        for j in range(n):
            W[i, j] = K[i, j] / d
    print(f' t={float(t):<8g} eps={float(e):.1e}  rho(P)-1 = {float(rho(W * (mp.eye(n) - t * Bc)) - 1):.3e}')
