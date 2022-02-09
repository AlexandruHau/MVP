import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import pandas as pd

# The number of lattices on one row / column. The overall lattice
# will have the dimensions N * N
N = 50
T = 2

# Declare using numpy an initial matrix lattice
spin_lattice = np.random.rand(50,50)

# All the numbers greater than 0.5 are set to 1 and all the elements
# smaller than 0.5 are converted to 0. This is done in order to get
# an initial Ising model of lattice with up and down spins
spin_lattice[spin_lattice > 0.5] = 1
spin_lattice[spin_lattice < 0.5] = -1

# Function for working out the overall energy of the lattice by analyzing
# the nearest-neighbours of each spin
def LatticeEnergy(S):
	ArrayEnergy = 0
	for i in range(N):
		for j in range(N):
			ArrayEnergy -= ( S[(i+1)%N][j] * S[i][j] + S[i][j] * S[i][(j+1)%N] )
	return ArrayEnergy
	
def LatticeMagnetization(S):
	Magnetization = 0
	for i in range(N):
		for j in range(N):
			Magnetization += S[i][j]
	return Magnetization

# Calculate the change in energy for the lattice by calculating only the spin
# product between the nearest neighbours	
def deltaE(S,i,j):
	dE = 2 * S[i][j] * (S[(i-1)%N][j] + S[(i+1)%N][j] + S[i][(j-1)%N] + S[i][(j+1)%N])
	return dE
	
# Function representing the Metropolis Algorithm implemented within the Glauber Algorithm.
# After one random state is flipped, check the probabilty of switch by applying the 
# Boltzmann weights
def SpinFlip(S, dE, T, i, j):
	threshold = np.exp(-dE/T)
	test_no = np.random.rand()
	if(test_no < threshold):
		S[i][j] = -S[i][j]
		# print("i row: " + str(i))
		# print("j column: " + str(j)) 
		
# Loop for Glauber Dynamics. Initialize a counter and perform this loop infinitely. We can
# only stop the loop by killing the program in the terminal through Ctrl+C command
def GlauberUpdate(T):
	dummy_variable = 0
	counter = 0
	
	# Calculate the initial energy of the Ising Model lattice
	E0 = LatticeEnergy(spin_lattice)
	
	# Initialize lists for the measurement and the overall magnetization. These two numpy
	# arrays will be plotted in the end on Matplotlib
	magnetizationLattice = []
	magnetizationSquared = []
	energy = []
	energySquared = []
	time = []
	while(counter < 1):
		for lx in range(N):
			for ly in range(N):
	
				# Select random coordinates for the cell
				i = np.random.randint(0,N)
				j = np.random.randint(0,N)
		
				# Work out the energy difference by flipping the spin in the lattice
				dE = deltaE(spin_lattice, i, j)
				
				dummy_variable += 1
		
				# Perform now the Metropolis Algorithm for spin flipping
				SpinFlip(spin_lattice, dE, T, i, j)
		
				# Use the counter variable for plotting the function
				# if(dummy_variable % 3000 == 0):
				# 	image = plt.imshow(spin_lattice, animated=True)
				# 	plt.draw()
				# 	plt.pause(0.01)
					
				# Function for plotting the magnetization of the lattice
				# as function of temperature
				if((dummy_variable % 25000 == 0) and (dummy_variable > 250000)):
					time.append(dummy_variable)
					M = LatticeMagnetization(spin_lattice)
					E= LatticeEnergy(spin_lattice)
					magnetizationLattice.append(np.abs(M))
					magnetizationSquared.append((np.abs(M))**2)
					energy.append(E)
					energySquared.append(E**2)
				
				# Exit the function
				if(dummy_variable > 25250000):
					return time, np.mean(magnetizationLattice), np.mean(magnetizationSquared), np.mean(energy), 						np.mean(energySquared)		
					counter = 1

# Now make a Python code for the Kawasaki Update step. In the Kawasaki Model, instead of
# spin flipping, two random spins should be flipped. Afterwards, the Metropolis Algorithm
# is implemented in order to check whether we accept or not the spin swap. 

# Introduce, first of all, a helper function for the spin swap in the scenario where the
# Metropolis Algorithm yields a successful result
def SpinSwap(S, dE, T, i1, j1, i2, j2):
	threshold = np.exp(-dE/T)
	test_no = np.random.rand()
	if(test_no < threshold):
		auxiliary_value = S[i1][j1]
		S[i1][j1] = S[i2][j2]
		S[i2][j2] = auxiliary_value

def KawasakiUpdate():
	dummy_variable = 0
	counter = 0
	while(counter < 1):
		for lx in range(N):
			for ly in range(N):
				
				# Select once again the random coordinates. NOTE: This time, 
				# we select four random coordinates for the two cells
				i1 = np.random.randint(0,N)
				j1 = np.random.randint(0,N)
				i2 = np.random.randint(0,N)
				j2 = np.random.randint(0,N)
				
				# Calculate the distance between the two points
				distance = np.sqrt((i1-i2)**2 + (j1-j2)**2)
				
				# Only analyse the case where the randomly selected spins are
				# opposite, otherwise there would be no change in the system.
				# Moreover, for the spin swap and energy change, take the first
				# scenario where the cells are not neighbours
				if((spin_lattice[i1][j1] != spin_lattice[i2][j2]) and (distance > 1)):	
					
					# Calculate the overall energy change from the swap 
					# of the spins
					dE = deltaE(spin_lattice, i1, j1) + deltaE(spin_lattice, i2, j2)
					
					# Perform now the Metropolis Algorithm for the spin swap
					SpinSwap(spin_lattice, dE, T, i1, j1, i2, j2)
				
				# Now analyse the case where the spins are swapped and the two cells
				# are neighbours. Double counting needs to be avoided for the neighboring
				# cells
				elif((spin_lattice[i1][j1] != spin_lattice[i2][j2]) and (distance == 1)):	
					
					# Calculate the overall energy change from the swap 
					# of the spins
					dE = deltaE(spin_lattice, i1, j1) + deltaE(spin_lattice, i2, j2) - 2
					
					# Perform now the Metropolis Algorithm for the spin swap
					SpinSwap(spin_lattice, dE, T, i1, j1, i2, j2)
					
				# Use the counter variable for plotting the function
				dummy_variable += 1
				if(dummy_variable % 1000 == 0):
					image = plt.imshow(spin_lattice, animated=True)
					plt.draw()
					plt.pause(0.01)

temperature = []
magnetization = []
susceptibility = []
energy = []
heatCapacity = []
count = 0
T = 1
while(T < 3):
	t1 = time.time()
	t, m, msq, E, Esq = GlauberUpdate(T)
	t2 = time.time()
	chi_M = (msq - m**2) / (N**2 * T)
	C_v = (Esq - E**2) / (N**2 * T**2)
	temperature.append(T)
	magnetization.append(m)
	energy.append(E)
	susceptibility.append(chi_M)
	heatCapacity.append(C_v)
	print("step " + str(count) + " finished with the taken time: " + str((t2-t1)/60) + " minutes")
	print("mean magnetization: " + str(m))
	print("mean lattice energy: " + str(E))
	print("magnetic susceptibility: " + str(chi_M) + "\n")
	count += 1
	T += 0.1

print("Temperature lattice: ")
print(temperature)
print("Magnetization: ")
print(magnetization)
print("Magnetic Susceptibility: ")
print(susceptibility)
print("Lattice Energy: ")
print(energy)
print("Heat Capacity: ")
print(heatCapacity)	

finalArray = np.asarray(temperature, magnetization, susceptibility, energy, heatCapacity)
pd.DataFrame(finalArray).to_csv("MCIsingModel.csv")
'''
fig, axis = plt.subplots(2)
print(len(temperature))
print(len(magnetization))
print(len(susceptibility))
print(len(energy))	
axis[0].scatter(temperature, magnetization)
axis[1].scatter(temperature, chi)
plt.show()
'''
