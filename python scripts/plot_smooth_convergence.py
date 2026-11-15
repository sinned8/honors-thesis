import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

RESOLUTIONS = (200, 400, 800, 1600)
RHO0 = 1.0e-3
P0 = 1.0e6


def load_solution(path, n):
    data = np.genfromtxt(path, delimiter=",", names=True)

    x = np.atleast_1d(data["x"])
    rho = np.atleast_1d(data["Density"])
    velocity = np.atleast_1d(data["Velocity"])
    pressure = np.atleast_1d(data["Pressure"])

    return x, rho, velocity, pressure


def average_to_coarse(values, n_coarse):
    n_fine = len(values)
    factor = n_fine // n_coarse
    return values.reshape(n_coarse, factor).mean(axis=1)


def print_convergence(variable, label, diff1, diff2, scale=None):
    e1 = np.mean(np.abs(diff1))
    e2 = np.mean(np.abs(diff2))
    linf1 = np.max(np.abs(diff1))
    linf2 = np.max(np.abs(diff2))

    print(f"\n{variable} self-convergence: {label}")
    print(f"L1 differences:   {e1:.8e}, {e2:.8e}")
    print(f"Linf differences: {linf1:.8e}, {linf2:.8e}")
    if scale is not None:
        print(f"L1 / reference scale: {e1 / scale:.8e}, {e2 / scale:.8e}")
    if e1 > 0 and e2 > 0:
        print(f"Observed L1 order: log2(E1/E2) = {np.log2(e1 / e2):.4f}")
    else:
        print("Observed L1 order is undefined (zero difference)")


def compare_variable(solutions, column, variable,
                     filename, output_dir, scale=None):
    # Triplet 1: 200, 400, 800 -- common N=200 grid
    x200 = solutions[200][0]
    q200 = solutions[200][column]
    q400_on_200 = average_to_coarse(solutions[400][column], 200)
    q800_on_200 = average_to_coarse(solutions[800][column], 200)

    diff_200_400 = q200 - q400_on_200
    diff_400_800_on_200 = q400_on_200 - q800_on_200

    print_convergence(variable, "200 - 400 - 800 (on N=200 grid)",
                      diff_200_400, diff_400_800_on_200, scale)

    # Triplet 2: 400, 800, 1600 -- common N=400 grid
    x400 = solutions[400][0]
    q400 = solutions[400][column]
    q800_on_400 = average_to_coarse(solutions[800][column], 400)
    q1600_on_400 = average_to_coarse(solutions[1600][column], 400)

    diff_400_800_on_400 = q400 - q800_on_400
    diff_800_1600 = q800_on_400 - q1600_on_400

    print_convergence(variable, "400 - 800 - 1600 (on N=400 grid)",
                      diff_400_800_on_400, diff_800_1600, scale)

    # Density is divided by background density. pressure and velocity use absolute differences.
    divisor = scale if scale is not None else 1.0
    ylabel = ("Absolute density difference / rho0" if variable == "Density"
              else f"Absolute {variable.lower()} difference")
    suffix = " / rho0" if variable == "Density" else ""

    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    axes[0].semilogy(x200, np.abs(diff_200_400) / divisor,
                     label="|N200 - N400|" + suffix)
    axes[0].semilogy(x200, np.abs(diff_400_800_on_200) / divisor,
                     label="|N400 - N800|" + suffix)
    axes[0].set(title="Triplet 1: comparison on N = 200 grid", ylabel=ylabel)

    axes[1].semilogy(x400, np.abs(diff_400_800_on_400) / divisor,
                     label="|N400 - N800|" + suffix)
    axes[1].semilogy(x400, np.abs(diff_800_1600) / divisor,
                     label="|N800 - N1600|" + suffix)
    axes[1].set(title="Triplet 2: comparison on N = 400 grid", xlabel="x", ylabel=ylabel)

    for ax in axes:
        ax.legend()
        ax.grid(alpha=0.25)

    fig.tight_layout()
    # plot_path = output_dir / filename
    # fig.savefig(plot_path, dpi=200)
    # print(f"Saved: {plot_path}")


def main():
    parser = argparse.ArgumentParser()
    default_output = Path(__file__).resolve().parent.parent / "outputs"
    parser.add_argument("--output-dir", type=Path, default=default_output)
    args = parser.parse_args()
    output_dir = args.output_dir.expanduser().resolve()

    missing = [output_dir / f"smooth_N{n}.csv" for n in RESOLUTIONS
               if not (output_dir / f"smooth_N{n}.csv").is_file()]
    if missing:
        print("Missing CSV files:")
        for path in missing:
            print("  ", path)
        return

    # each entry is (x, density, velocity, pressure)
    solutions = {}
    for n in RESOLUTIONS:
        solutions[n] = load_solution(output_dir / f"smooth_N{n}.csv", n)

    #
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for n in RESOLUTIONS:
        x, _, _, pressure = solutions[n]
        ax.plot(x, pressure - P0, label=f"N = {n}", linewidth=1.4)
    ax.set(xlabel="x", ylabel="Pressure perturbation p - p0",
           title="Smooth-flow test: pressure perturbation")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    # pressure_plot = output_dir / "smooth_pressure_comparison.png"
    # fig.savefig(pressure_plot, dpi=200)
    # print(f"Saved: {pressure_plot}")

    compare_variable(solutions, 1, "Density", "smooth_density_self_convergence.png",
                     output_dir, scale=RHO0)
    compare_variable(solutions, 2, "Velocity", "smooth_velocity_self_convergence.png",
                     output_dir)
    compare_variable(solutions, 3, "Pressure", "smooth_pressure_self_convergence.png",
                     output_dir)

    plt.show()


if __name__ == "__main__":
    main()
