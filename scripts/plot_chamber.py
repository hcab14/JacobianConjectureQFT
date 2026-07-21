"""Figure: the C*-reduced chamber geometry of the counterexample.

Draws, in the invariant plane (u, w) = (a c^2, b c):
  * the wall {P2 = 0} (cuspidal cubic, the reduced non-properness locus),
  * the N=3 chamber {P2 < 0} (where the real map is 3-to-1),
  * the cusp (4/27, 4/3) = the C*-orbit missed by F (empty fiber),
  * the cuspidal tangent = the reduced x-collision line {D0 = 0},
  * the image of the tree-expansion ray J = t(1,2,3), whose first wall
    crossing at |t| ~ 0.302 is the convergence radius.

Output: docs/img/chamber_geometry.png
"""

import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

U0, W0 = 4 / 27, 4 / 3

fig, (ax, ax2) = plt.subplots(
    1, 2, figsize=(12.5, 6.2), gridspec_kw={"width_ratios": [1.15, 1]})

for axis, (ulim, wlim), title in (
    (ax, ((-4.0, 4.0), (-6.0, 2.4)), "the chamber and its wall"),
    (ax2, ((-0.32, 0.42), (0.35, 1.75)), "zoom: the cuspidal horn"),
):
    uu = np.linspace(*ulim, 900)
    ww = np.linspace(*wlim, 900)
    GU, GW = np.meshgrid(uu, ww)
    P2 = 27 * GU**2 + 16 * GU - 18 * GU * GW + GW**3 - GW**2

    axis.contourf(GU, GW, P2, levels=[-1e9, 0], colors=["#aecbe8"], alpha=0.85)
    axis.contour(GU, GW, P2, levels=[0], colors="#1f4e79", linewidths=1.8)

    # cuspidal tangent = reduced D0-line: 27u - 9w + 8 = 0
    axis.plot(uu, 3 * uu + 8 / 9, color="#c05020", lw=1.2, ls="--",
              label=r"$\{D_0=0\}$: cuspidal tangent" if axis is ax else None)

    axis.plot([U0], [W0], marker="o", ms=7, mfc="#d62728", mec="k", zorder=5)
    axis.annotate("cusp $(4/27,\\,4/3)$\n(empty fiber: all 3 sheets\nat infinity)",
                  (U0, W0), textcoords="offset points", xytext=(12, 8),
                  fontsize=9)

    # image of the ray J = t(1,2,3):  (u,w) = (9 t^3, 6 t^2)
    tt = np.linspace(-0.75, 0.75, 601)
    axis.plot(9 * tt**3, 6 * tt**2, color="#2ca02c", lw=1.6,
              label=r"ray $J=t(1,2,3)$: $(9t^3,\,6t^2)$" if axis is ax else None)
    ts = -0.302028  # nearest wall crossing of the ray (t negative)
    axis.plot([9 * ts**3], [6 * ts**2], marker="s", ms=6, mfc="#2ca02c",
              mec="k", zorder=5)
    if axis is ax2:
        axis.annotate("first wall crossing $t^*\\approx -0.302$:\n"
                      "$|t^*|$ = radius of convergence\nof the tree series",
                      (9 * ts**3, 6 * ts**2), textcoords="offset points",
                      xytext=(-4, -52), fontsize=9)

    axis.set_xlim(*ulim)
    axis.set_ylim(*wlim)
    axis.set_xlabel(r"$u = a\,c^2$")
    axis.set_ylabel(r"$w = b\,c$")
    axis.set_title(title, fontsize=11)
    axis.axhline(0, color="0.6", lw=0.5)
    axis.axvline(0, color="0.6", lw=0.5)

ax.text(-3.7, -5.5,
        r"$N=3$ chamber: $P_2<0$" + "\n" + r"($p<0$, three real preimages)",
        fontsize=10, color="#1f4e79")
ax.text(2.3, 0.4, r"$N=1$: $P_2>0$", fontsize=10, color="0.25")
ax.legend(loc="lower right", fontsize=9, framealpha=0.9)

fig.suptitle(
    r"C$^*$-reduced source space of the Alpöge–Mathew map: "
    r"the wall $\{p=0\}$ is a cuspidal cubic $P_2(u,w)=27u^2+16u-18uw+w^3-w^2$",
    fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.95))

out = pathlib.Path(__file__).resolve().parent.parent / "docs" / "img"
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "chamber_geometry.png", dpi=160)
print(f"wrote {out / 'chamber_geometry.png'}")
