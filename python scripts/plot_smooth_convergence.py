
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

RESOLUTIONS = (200, 400, 800, 1600, 3200)
TRIPLETS = ((200, 400, 800), (400, 800, 1600), (800, 1600, 3200))
RHO0 = 1.0e-3
P0 = 1.0e6


def load_solution(path, n):
    data = np.genfromtxt(path, delimiter=",", names=True)
    columns = ("x", "Density", "Velocity", "Pressure")
    if data.dtype.names is None or any(name not in data.dtype.names for name in columns):
        raise ValueError(f"{path.name}: expected columns x,Density,Velocity,Pressure")

    x = np.atleast_1d(data["x"])
    rho = np.atleast_1d(data["Density"])
    velocity = np.atleast_1d(data["Velocity"])
    pressure = np.atleast_1d(data["Pressure"])

    if len(x) != n:
        raise ValueError(f"{path.name}: expected {n} rows, found {len(x)}")
    if not all(np.isfinite(a).all() for a in (x, rho, velocity, pressure)):
        raise ValueError(f"{path.name}: NaN or infinite value found")
    if not np.allclose(x, (np.arange(n) + 0.5)/n, rtol=0, atol=1e-10):
        raise ValueError(f"{path.name}: unexpected x cell centers")
    return x, rho, velocity, pressure


def average_to_coarse(values, n_coarse):
    n_fine = len(values)
    if n_fine % n_coarse != 0:
        raise ValueError("Fine resolution must be a multiple of coarse resolution")
    return values.reshape(n_coarse, n_fine // n_coarse).mean(axis=1)


def compare_variable(solutions, column, variable, filename, output_dir, scale=None):
    fig, axes = plt.subplots(3, 1, figsize=(9, 11), sharex=True)
    divisor = scale if scale is not None else 1.0
    ylabel = "Absolute density difference / rho0" if variable == "Density" else f"Absolute {variable.lower()} difference"
    suffix = " / rho0" if variable == "Density" else ""

    for ax, (coarse, middle, fine) in zip(axes, TRIPLETS):
        x = solutions[coarse][0]
        q_coarse = solutions[coarse][column]
        q_middle = average_to_coarse(solutions[middle][column], coarse)
        q_fine = average_to_coarse(solutions[fine][column], coarse)

        diff1 = q_coarse - q_middle
        diff2 = q_middle - q_fine
        e1 = np.mean(np.abs(diff1))
        e2 = np.mean(np.abs(diff2))
        linf1 = np.max(np.abs(diff1))
        linf2 = np.max(np.abs(diff2))

        print(f"\n{variable} self-convergence: {coarse} - {middle} - {fine} (on N={coarse} grid)")
        print(f"L1 differences:   {e1:.8e}, {e2:.8e}")
        print(f"Linf differences: {linf1:.8e}, {linf2:.8e}")
        if scale is not None:
            print(f"L1 / reference scale: {e1/scale:.8e}, {e2/scale:.8e}")
        if e1 > 0 and e2 > 0:
            print(f"Observed L1 order: log2(E1/E2) = {np.log2(e1/e2):.4f}")
            print(f"Observed Linf order: log2(E1/E2) = {np.log2(linf1/linf2):.4f}" if linf2 > 0 else "Observed Linf order undefined")
        else:
            print("Observed order undefined (zero difference)")

        # The curves show pointwise differences; L1/Linf above are computed before scaling.
        ax.semilogy(x, np.maximum(np.abs(diff1) / divisor, 1e-18),
                    label=f"|N{coarse} - N{middle}|" + suffix)
        ax.semilogy(x, np.maximum(np.abs(diff2) / divisor, 1e-18),
                    label=f"|N{middle} - N{fine}|" + suffix)
        ax.set(title=f"Comparison on N = {coarse} grid", ylabel=ylabel)
        ax.legend()
        ax.grid(alpha=0.25)

    axes[-1].set_xlabel("x")
    fig.tight_layout()
    # path = output_dir / filename
    # fig.savefig(path, dpi=200)
    # print(f"Saved: {path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).resolve().parent.parent / "outputs")
    args = parser.parse_args()
    output_dir = args.output_dir.expanduser().resolve()

    solutions = {}
    for n in RESOLUTIONS:
        path = output_dir / f"smooth_N{n}.csv"
        if not path.is_file():
            print(f"Missing: {path}")
            return
        solutions[n] = load_solution(path, n)

    fig, ax = plt.subplots(figsize=(9, 4.5))
    for n in RESOLUTIONS:
        x, _, _, pressure = solutions[n]
        ax.plot(x, pressure - P0, label=f"N = {n}", linewidth=1.3)
    ax.set(xlabel="x", ylabel="Pressure perturbation p - p0",
           title="Smooth-flow test: pressure perturbation")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    # fig.savefig(output_dir / "smooth_pressure_comparison.png", dpi=200)

    compare_variable(solutions, 1, "Density", "smooth_density_self_convergence.png",
                     output_dir, scale=RHO0)
    compare_variable(solutions, 2, "Velocity", "smooth_velocity_self_convergence.png",
                     output_dir)
    compare_variable(solutions, 3, "Pressure", "smooth_pressure_self_convergence.png",
                     output_dir)
    plt.show()


if __name__ == "__main__":
    main()
