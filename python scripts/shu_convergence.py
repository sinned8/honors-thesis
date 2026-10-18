

from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs"
REF = ROOT / "external" / "SHU_RHO_PROFILE_EXACT.mat"
LEVELS = [200, 400, 800, 1600]
X_MIN, X_MAX = -5.0, 5.0

#This is the a file to test the convergence of various error norms against the supplied reference
#solution from https://github.com/sulli72/1D_EULER/tree/master


def reference_cell_averages(x_ref, rho_ref, edges):
    averages = np.empty(len(edges) - 1)

    for i in range(len(averages)):
        left, right = edges[i], edges[i + 1]
        # Include the reference points strictly inside this cell.
        inside = (x_ref > left) & (x_ref < right)
        x = np.r_[left, x_ref[inside], right]
        rho = np.interp(x, x_ref, rho_ref)
        integral = np.sum(0.5 * (rho[:-1] + rho[1:]) * np.diff(x))
        averages[i] = integral / (right - left)

    return averages


def load_numerical_csv(filename, n):
    data = np.genfromtxt(filename, delimiter=",", names=True)
    x = np.atleast_1d(data["x"])
    rho = np.atleast_1d(data["Density"])

    if len(x) != n or not np.all(np.isfinite(rho)):
        raise ValueError(f"{filename}: expected {n} finite density values")

    dx = (X_MAX - X_MIN) / n
    expected_x = X_MIN + (np.arange(n) + 0.5) * dx
    if not np.allclose(x, expected_x, rtol=0.0, atol=1e-8):
        raise ValueError(f"{filename}: x positions don't match {n} FV cell centers")

    return x, rho


def main():
    ref = loadmat(REF)
    x_ref = np.asarray(ref["xv"], dtype=float).ravel()
    rho_ref = np.asarray(ref["r"], dtype=float).ravel()

    order = np.argsort(x_ref)
    x_ref, rho_ref = x_ref[order], rho_ref[order]

    if len(x_ref) != len(rho_ref) or not np.all(np.diff(x_ref) > 0):
        raise ValueError("Reference coordinates/density are invalid")
    if not (np.all(np.isfinite(x_ref)) and np.all(np.isfinite(rho_ref))):
        raise ValueError("Non-finite values found in the reference")
    if x_ref[0] > X_MIN + 1e-8 or x_ref[-1] < X_MAX - 1e-8:
        raise ValueError("The reference must cover the full [-5, 5] domain")

    results = []
    previous = None

    for n in LEVELS:
        filename = OUT / f"shu_N{n}.csv"
        x, rho = load_numerical_csv(filename, n)
        edges = np.linspace(X_MIN, X_MAX, n + 1)
        rho_ref_avg = reference_cell_averages(x_ref, rho_ref, edges)
        error = np.abs(rho - rho_ref_avg)

        l1 = np.mean(error)
        l2 = np.sqrt(np.mean(error ** 2))
        linf = np.max(error)
        errors = [l1, l2, linf]

        if previous is None:
            rates = [np.nan, np.nan, np.nan]
        else:
            rates = [np.log(old / new) / np.log(n / last_n)
                     for old, new in zip(previous, errors)]

        results.append([n, *errors, *rates])
        previous, last_n = errors, n

        # Keep an aligned numerical/reference profile for later inspection.
        np.savetxt(
            OUT / f"shu_aligned_N{n}.csv",
            np.column_stack((x, rho, rho_ref_avg, error)),
            delimiter=",", comments="",
            header="x,Density_Numerical,Density_ReferenceAvg,Absolute_Error",
            fmt="%.15g"
        )

    result_file = OUT / "shu_convergence.csv"
    with open(result_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["N", "L1", "L2", "Linf", "L1_Order", "L2_Order", "Linf_Order"])
        writer.writerows(results)

    print("\nShu-Osher density errors (WENO3 + SSP-RK2 vs WENO5 reference)")
    print(f"{'N':>6} {'L1':>12} {'L2':>12} {'Linf':>12} {'p_L1':>8} {'p_L2':>8} {'p_Linf':>8}")
    for n, l1, l2, linf, p1, p2, pinf in results:
        print(f"{n:6d} {l1:12.5e} {l2:12.5e} {linf:12.5e} "
              f"{p1:8.3f} {p2:8.3f} {pinf:8.3f}")

    ns = [row[0] for row in results]
    plt.figure(figsize=(8, 5))
    for index, label in [(1, "L1"), (2, "L2"), (3, "Linf")]:
        plt.loglog(ns, [row[index] for row in results], "o-", label=label)
    plt.xlabel("Number of cells (N)")
    plt.ylabel("Density error relative to WENO5 reference")
    plt.title("Shu-Osher density convergence, t = 1.8")
    plt.xticks(ns, [str(n) for n in ns])
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / "shu_error_vs_resolution.png", dpi=200)
    plt.close()

    # High-resolution overlay focused on the shock/oscillation region.
    fine = np.genfromtxt(OUT / "shu_N1600.csv", delimiter=",", names=True)
    plt.figure(figsize=(10, 5))
    plt.plot(x_ref, rho_ref, label="WENO5 reference (4000 points)")
    plt.plot(fine["x"], fine["Density"], "--", label="WENO3 + SSP-RK2 (N=1600)")
    plt.xlim(0.5, 2.5)
    plt.xlabel("x")
    plt.ylabel("Density")
    plt.title("Shu-Osher shock-region comparison, t = 1.8")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / "shu_reference_overlay_zoom.png", dpi=200)
    plt.close()

    print(f"\nSaved table: {result_file}")
    print("Saved plots: outputs/shu_error_vs_resolution.png and outputs/shu_reference_overlay_zoom.png")
    print("Saved aligned cell data: outputs/shu_aligned_N*.csv")


if __name__ == "__main__":
    main()
