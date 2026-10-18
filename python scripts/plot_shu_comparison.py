import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat

# WENO3 + SSP-RK2 output
data = pd.read_csv("outputs/shu_N1600.csv")

x = data["x"].to_numpy()
rho = data["Density"].to_numpy()

# WENO5 high-resolution reference thanks to https://github.com/sulli72/1D_EULER/tree/master
reference = loadmat("external/SHU_RHO_PROFILE_EXACT.mat")

x_ref = reference["xv"].flatten()
rho_ref = reference["r"].flatten()

#plotting both density profiles
plt.figure(figsize=(11, 5))

plt.plot(x_ref, rho_ref, "k-",
         label="WENO5 Reference (N=4000)")
plt.plot(x, rho, "r--",
         label="WENO3 + SSP-RK2 (N=1600)")

plt.xlabel("x")
plt.ylabel("Density")
plt.title("Shu-Osher Density Comparison, t = 1.8")

plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

plt.show()