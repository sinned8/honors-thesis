import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np

from sod_exact import exact_sod_solution,calculate_errors


x_exact = np.linspace(0.0, 1.0, 1000)

t = 0.2

#calculating exact solution
rho_exact, u_exact, p_exact = exact_sod_solution(
    x_exact,
    t,
    x0=0.5,
    gamma=1.4
)




csv_file = sys.argv[1]
data = pd.read_csv(csv_file)

x_num = data["x"].to_numpy()
rho_num = data["Density"].to_numpy()
u_num = data["Velocity"].to_numpy()
p_num = data["Pressure"].to_numpy()


N = len(x_num)
dx = 1.0 / N
#exact density at cell centers
rho_exact_cells, _, _ = exact_sod_solution(
    x_num,
    t,
    x0=0.5,
    gamma=1.4
)
#calculating and displaying density error
L1, L2, Linf = calculate_errors(rho_num, rho_exact_cells, dx)
print(f"N = {N}")
print(f"L1 = {L1:.8e}")
print(f"L2 = {L2:.8e}")
print(f"Linf = {Linf:.8e}")

fig, axes = plt.subplots(3, 1, figsize=(8, 9))

# Density
axes[0].plot(x_exact, rho_exact, label="Exact")
axes[0].plot(x_num, rho_num, "o", markersize=2, label="WENO3 + Rusanov + SSP-RK2")
axes[0].set_ylabel("Density")
axes[0].legend(loc="lower left")

# Zoom in on the contact discontinuity to see small differences
zoom = axes[0].inset_axes([0.60, 0.53, 0.37, 0.39])
zoom.plot(x_exact, rho_exact)
zoom.plot(x_num, rho_num, "o", markersize=2)
zoom.set_xlim(0.65, 0.71)
zoom.set_ylim(0.23, 0.46)
zoom.set_title("Contact discontinuity", fontsize=9)
zoom.tick_params(labelsize=8)
axes[0].indicate_inset_zoom(zoom, edgecolor="gray")

# Velocity
axes[1].plot(x_exact, u_exact, label="Exact")
axes[1].plot(x_num, u_num, "o", markersize=2, label="WENO3 + Rusanov + SSP-RK2")
axes[1].set_ylabel("Velocity")
axes[1].legend()

# Pressure
axes[2].plot(x_exact, p_exact, label="Exact")
axes[2].plot(x_num, p_num, "o", markersize=2, label="WENO3 + Rusanov + SSP-RK2")
axes[2].set_ylabel("Pressure")
axes[2].set_xlabel("x")
axes[2].legend()

fig.suptitle(f"Sod Problem — WENO3 + Rusanov + SSP-RK2, N = {N}, t = {t}")

plt.tight_layout()
plt.show()

