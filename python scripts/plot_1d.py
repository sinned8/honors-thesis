import pandas as pd
import sys
import matplotlib.pyplot as plt


csv_file = sys.argv[1]
df = pd.read_csv(csv_file)

plt.figure(1)
plt.plot(df["x"], df["Density"])
plt.xlabel("x")
plt.ylabel("Density")
plt.title("Sod Shock Tube - Density at t = 0.2")
plt.grid()

plt.figure(2)
plt.plot(df["x"], df["Velocity"])
plt.xlabel("x")
plt.ylabel("Velocity")
plt.title("Sod Shock Tube - Velocity at t = 0.2")
plt.grid()

plt.figure(3)
plt.plot(df["x"], df["Pressure"])
plt.xlabel("x")
plt.ylabel("Pressure")
plt.title("Sod Shock Tube - Pressure at t = 0.2")
plt.grid()


plt.show()