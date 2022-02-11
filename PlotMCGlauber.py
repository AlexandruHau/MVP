# Python file for the plot of the data extracted from
# the Glauber Monte Carlo Dynamics
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Read the data from the generated simulations using the
# lattice Ising Model
df = pd.read_csv("MCIsingModel.csv")

# Take all the arrays in numpy forms
temperature = df["Temperature"].to_numpy()
magnetization = df["Magnetization"].to_numpy()
energy = df["Lattice Energy"].to_numpy()
cV = df["Heat Capacity"].to_numpy()
chi = df["Magnetic Susceptibility"].to_numpy()
err_cV = df["Error on heat capacity"].to_numpy()
err_chi = df["ES"].to_numpy()

# Plot commands on matplotlib
# All the four temperature evolution are shown on the
# same mother plot
fig, axis = plt.subplots(2,2)
fig.suptitle("Plots on the Ising Lattice - Glauber Dynamics")
axis[0,0].scatter(temperature, magnetization)
axis[0,0].set_xlabel("Temperature")
axis[0,0].set_ylabel("Magnetization lattice")
axis[0,1].scatter(temperature, energy)
axis[0,1].set_xlabel("Temperature")
axis[0,1].set_ylabel("Energy lattice")
axis[1,0].errorbar(temperature, chi, yerr=err_chi, fmt="o")
axis[1,0].set_xlabel("Temperature")
axis[1,0].set_ylabel("Magnetic susceptibility")
axis[1,1].errorbar(temperature, cV, yerr=err_cV, fmt="o")
axis[1,1].set_xlabel("Temperature")
axis[1,1].set_ylabel("Heat capacity")
plt.tight_layout()
plt.show()
