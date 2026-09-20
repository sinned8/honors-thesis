import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np

from sod_exact import exact_sod_solution


x_exact = np.linspace(0.0, 1.0, 1000)

t = 0.2

rho_exact, u_exact, p_exact = exact_sod_solution(
    x_exact,
    t,
    x0=0.5,
    gamma=1.4
)




csv_file = sys.argv[1]
data = pd.read_csv(csv_file)

x_num = data["x"]
rho_num = data["Density"]
u_num = data["Velocity"]
p_num = data["Pressure"]

fig, axes = plt.subplots(3, 1, figsize=(8, 9))

# Density
axes[0].plot(x_exact, rho_exact, label="Exact")
axes[0].plot(x_num, rho_num, "o", markersize=2, label="Rusanov")
axes[0].set_ylabel("Density")
axes[0].legend()

# Velocity
axes[1].plot(x_exact, u_exact, label="Exact")
axes[1].plot(x_num, u_num, "o", markersize=2, label="Rusanov")
axes[1].set_ylabel("Velocity")
axes[1].legend()

# Pressure
axes[2].plot(x_exact, p_exact, label="Exact")
axes[2].plot(x_num, p_num, "o", markersize=2, label="Rusanov")
axes[2].set_ylabel("Pressure")
axes[2].set_xlabel("x")
axes[2].legend()

plt.tight_layout()
plt.show()

