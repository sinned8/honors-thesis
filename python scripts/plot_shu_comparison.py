import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat
from pathlib import Path




##**********************************************************|
#           Plots FV Shu Osher reference results against    |
#           various spatial resolutions of my solver        |
##**********************************************************|



root = Path(__file__).resolve().parent.parent
output_dir = root / "outputs"

# WENO5 high-resolution reference thanks to https://github.com/sulli72/1D_EULER/tree/master
reference = loadmat(root / "external" / "SHU_RHO_PROFILE_EXACT.mat")
x_ref = reference["xv"].flatten()
rho_ref = reference["r"].flatten()
resolutions = [200,400,800,1600, 4000, 6400]

# Load all curves first
curves = {}
for n in resolutions:
    f = output_dir / f"shu_N{n}.csv"
    if not f.exists():
        print(f"Missing: {f}")
        continue
    try:
        d = np.loadtxt(f, delimiter=",")
    except ValueError:
        d = np.loadtxt(f, delimiter=",", skiprows=1)
    curves[n] = (d[:, 0], d[:, 1])

colors = plt.cm.viridis(np.linspace(0.0, 0.8, len(curves)))  # dark -> light = coarse -> fine

fig, (ax_full, ax_zoom) = plt.subplots(2, 1, figsize=(11, 8.5))
axins = ax_zoom.inset_axes([0.03, 0.06, 0.30, 0.42])  # empty lower left of zoom panel used for added visibility

for ax in (ax_full, ax_zoom, axins):
    # Reference as a slim dashed gray halo underneath, so converged curves ride on top of it
    ax.plot(x_ref, rho_ref, color="0.55", linestyle="--", lw=2.5, alpha=0.55,
            label="Reference (WENO5)", zorder=1)
    for (n, (x, rho)), c in zip(curves.items(), colors):
        ax.plot(x, rho, color=c, lw=1.2, label=f"$N={n}$", zorder=2)
    ax.grid(alpha=0.3)

ax_full.axvspan(0.5, 2.5, color="0.9", zorder=0)   # shows what the zoom covers
ax_zoom.set_xlim(0.5, 2.5)
axins.set_xlim(1.75, 1.95)
axins.set_ylim(4.2, 4.75)
ax_zoom.indicate_inset_zoom(axins, edgecolor="0.4")

ax_full.set_title(r"Shu-Osher density, $t = 1.8$: full domain")
ax_zoom.set_title("Zoomed density oscillations")
ax_zoom.set_xlabel("$x$")
for ax in (ax_full, ax_zoom):
    ax.set_ylabel(r"$\rho$")

ax_full.legend(loc="lower left", ncol=2)  # empty region of the top panel

fig.tight_layout()
fig.savefig(root / "outputs" / "shu_osher_resolution_comparison.png", dpi=250)
plt.show()








