# Python file for the plot of the data extracted from
# the Kawasaki Monte Carlo Dynamics
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Read the data from the generated simulations using the
# lattice Ising Model
df = pd.read_csv("MCIsingModelKawasaki.csv")

# Take all the arrays in numpy forms
temperature = df["Temperature"].to_numpy()
energy = df["Lattice Energy"].to_numpy()
cV = df["Heat Capacity"].to_numpy()
err_cV = df["Error on heat capacity"].to_numpy()

# Plot commands on matplotlib. The two plots for lattice energy and
# heat capacity are shown on the same plot
fig, axis = plt.subplots(1,2)
fig.suptitle("Plots on the Ising Lattice - Kawasaki Dynamics")
axis[0].scatter(temperature, energy)
axis[0].set_xlabel("Temperature")
axis[0].set_ylabel("Magnetization lattice")
axis[1].errorbar(temperature, cV, yerr=err_cV, fmt="o")
axis[1].set_xlabel("Temperature")
axis[1].set_ylabel("Heat capacity")
plt.tight_layout()
plt.show()
