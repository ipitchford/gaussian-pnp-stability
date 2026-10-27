#!/usr/bin/env python3
"""Figure 1 of the revised manuscript: stability regions and growth rates.

Illustration only; every boundary drawn here is an exact formula proved in the
text (Delta of eq. (43), F of eq. (47)). Pinned plotting dependencies are
listed in README.md.
"""
from pathlib import Path
import os
import tempfile
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "pnp-matplotlib"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
INK, PLUM, GOLD, GREY = "#28212A", "#73365F", "#B47721", "#7d7378"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9.5,
                     "axes.spines.top": False, "axes.spines.right": False})


def delta(e, t):
    """Exact Delta(eps,t) of eq. (43): Delta < 0 <=> unstable (Theorem 4)."""
    return (e * (31/120 + 217*e/200)
            + t * (-133/3750000 + 1782779*e/10000000 + 969583*e*e/50000000)
            + t*t * (751389/12500000000 - 15987*e/1000000000 - 1582713*e*e/125000000000))


def F3(e, t):
    """Exact F_eps(t) of eq. (47) for the three-state family (Theorem 5)."""
    tau = 1 - e
    return e*(1+e)/2 + tau*(60*e-1)/240*t + 29*tau**2/7200*t*t


def first_loss(e):
    c = e*(31/120 + 217*e/200)
    b = -133/3750000 + 1782779*e/10000000 + 969583*e*e/50000000
    a = 751389/12500000000 - 15987*e/1000000000 - 1582713*e*e/125000000000
    return 2*c/(-b + np.sqrt(b*b - 4*a*c))


v = np.array([-8, -8, 17, -1.0])
A = (20*np.eye(4) + 245*np.ones((4, 4)) + np.outer(v, v))/1000
B = A @ A
K0 = np.array([[1, .9, .1, 0], [.9, 1, .1, 0], [.1, .1, 1, 0], [0, 0, 0, 1]])
H = np.array([[0, 1.5, .3, .2], [1.5, 0, .3, .2], [.3, .3, 0, .6], [.2, .2, .6, 0]])

fig, axes = plt.subplots(2, 2, figsize=(11, 8.4))
fig.subplots_adjust(left=.08, right=.98, bottom=.07, top=.93, wspace=.27, hspace=.42)

# (a) four-state stability region
ax = axes[0, 0]
E, T = np.meshgrid(np.geomspace(1e-9, 3e-5, 500), np.linspace(1e-4, 2, 500))
ax.contourf(E, T, (delta(E, T) < 0).astype(float), levels=[.5, 1.5], colors=[PLUM], alpha=.35)
ax.contour(E, T, delta(E, T), levels=[0], colors=[PLUM], linewidths=1.4)
ts = np.linspace(1e-4, 2, 400)
ax.plot(ts/1e5, ts, color=INK, lw=1, ls="--", label=r"family $\epsilon=t/10^5$ (Thm 1b)")
ax.axvline(1.6946633e-5, color=GOLD, lw=1, label=r"$\epsilon_{\rm crit}\approx1.69\times10^{-5}$")
ax.scatter([0.5312853/1e5], [0.5312853], color=INK, s=18, zorder=5)
ax.annotate(r"$t^*\approx0.5313$", (0.5312853/1e5, 0.5312853), (2.2e-8, 1.15), fontsize=9,
            arrowprops={"arrowstyle": "-", "color": INK})
ax.set_xscale("log"); ax.set(xlabel=r"weak-link strength $\epsilon$", ylabel="step $t$", ylim=(0, 2))
ax.set_title("(a) Four states: unstable region (shaded)", loc="left")
ax.legend(frameon=False, fontsize=8, loc="upper left")

# (b) growth per step along the family
ax = axes[0, 1]
steps = np.linspace(.002, .6, 300)
gp, gr = [], []
for t in steps:
    K = K0 + (t/1e5)*H
    W = K/K.sum(1)[:, None]
    P = W @ (np.eye(4) - t*B)
    R = np.eye(4) - W + np.linalg.solve(np.eye(4) + t*B, 2*W - np.eye(4))
    gp.append(np.abs(np.linalg.eigvals(P)).max() - 1)
    gr.append(np.abs(np.linalg.eigvals(R)).max() - 1)
ax.axhline(0, color=INK, lw=.7)
ax.plot(steps, 1e5*np.array(gp), color=PLUM, lw=2.4, label=r"$P(t)$ (PnP-ISTA)")
ax.plot(steps, 1e5*np.array(gr), color=GOLD, lw=1.4, ls=(0, (2, 2)), label=r"$R(t)$ (PnP-ADMM)")
ax.axvline(.5312853, color=GREY, lw=.8, ls=":")
ax.text(.535, 1.2, r"$t^*$", color=GREY)
ax.set(xlabel="step $t$", ylabel=r"$10^5\,[\rho-1]$ per iteration", xlim=(0, .6))
ax.set_title(r"(b) Growth per step is at most $1.7\times10^{-5}$", loc="left")
ax.legend(frameon=False, fontsize=8, loc="lower left")

# (c) first loss versus mixing gap
ax = axes[1, 0]
es = np.geomspace(1e-10, 1.6e-5, 300)
Q0 = np.array([[19/20, 1/20], [1/6, 5/6]]); Q1 = np.array([[7/20, -9/20], [1/10, -3/5]])
gs = np.array([1 - np.sort(np.linalg.eigvals((Q0 + e*Q1)/(1 + e)).real)[-1] for e in es])
ax.loglog(gs, first_loss(es), color=PLUM, lw=2.4, label="true first loss of stability")
ax.loglog(gs, gs/np.sqrt(2), color=GOLD, lw=1.8, label=r"guaranteed bound $g/\sqrt{\kappa_\pi}$ (Thm 2)")
ax.set(xlabel="mixing gap $g$", ylabel="step $t$")
ax.text(3e-9, 3e-6, r"ratio $\approx 8{,}600$", color=INK)
ax.set_title("(c) The guarantee is linear in $g$ but conservative", loc="left")
ax.legend(frameon=False, fontsize=8, loc="upper left")

# (d) three-state reentrant region (non-unit diagonal)
ax = axes[1, 1]
E3, T3 = np.meshgrid(np.linspace(1e-6, 2.2e-3, 400), np.linspace(1e-4, 2, 400))
ax.contourf(1e3*E3, T3, (F3(E3, T3) < 0).astype(float), levels=[.5, 1.5], colors=[GOLD], alpha=.4)
ax.contour(1e3*E3, T3, F3(E3, T3), levels=[0], colors=[GOLD], linewidths=1.4)
ax.axhline(.25, color=GREY, lw=.6, ls=":")
ax.set(xlabel=r"$10^3\,\epsilon$ (three-state family)", ylabel="step $t$", ylim=(0, 2))
ax.set_title("(d) Without unit diagonals: failure in $n=3$", loc="left")
ax.text(1.2, 1.5, "stable", color=INK); ax.text(.15, .45, "unstable", color=INK)

fig.savefig(HERE/"gaussian_stability_diagram.png", dpi=180)
fig.savefig(HERE/"gaussian_stability_diagram.pdf")
print("Figure saved: gaussian_stability_diagram.png/.pdf")
