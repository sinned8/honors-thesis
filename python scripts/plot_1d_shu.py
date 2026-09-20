import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np





csv_file = sys.argv[1]
N = sys.argv[2]
T_FINAL = sys.argv[3]


df = pd.read_csv(csv_file)

x = df["x"]
rho = df["Density"]
velocity = df["Velocity"]
pressure = df["Pressure"]



#debugging plot
fig, axes = plt.subplots(
    3, 1,
    figsize=(11, 8),
    sharex=True
)

axes[0].plot(x, rho, linewidth=1.2)
axes[0].set_ylabel("Density")
axes[0].grid(alpha=0.25)

axes[1].plot(x, velocity, linewidth=1.2)
axes[1].set_ylabel("Velocity")
axes[1].grid(alpha=0.25)

axes[2].plot(x, pressure, linewidth=1.2)
axes[2].set_ylabel("Pressure")
axes[2].set_xlabel("x")
axes[2].grid(alpha=0.25)

axes[2].set_xlim(-5, 5)

fig.suptitle(f"Shu–Osher Problem at {N}")

plt.tight_layout()
plt.show()




#comparison vs reference image plot
plt.figure(figsize=(10, 4.5))

plt.plot(
    x,
    rho,
    linewidth=1.3,
    label="Rusanov"
)

plt.xlabel("x")
plt.ylabel("Density")
plt.title(f"Shu–Osher Density at t = {T_FINAL}")
plt.xlim(-5, 5)

plt.grid(alpha=0.25)
plt.legend()

plt.tight_layout()
plt.show()