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

# The number of particles is given by the square of the number of
# cells per row
no_particles = N**2

# Declare using numpy an initial matrix lattice
spin_lattice = np.random.rand(50,50)

# Declare in advance the number of sweeps required before beginning
# the measurements as well as the number of sweeps required between
# measurements. One sweep consists of an iteration through the whole
# lattice
no_initial_sweeps = 100
no_sweeps_between = 10

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

# Function for working out the overall magnetization of the lattice by 
# adding up the spin of each particle	
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
		
def Bootstrap_Sampling(arr):
	# Bootstrap method for heat capacity in the Glauber Dynamics model
	# Declare firstly the number of samplings we would like to take.
	# Also declare the initial empty array of k samplings which will
	# be updated at each step
	# Declare the number of iterations for the bootstrap sampling
	k = 100
	T = 1

	# Declare the number of elements n the user wants to use for
	# each iteration
	n = len(arr)

	# Declare the number of particles
	no_particles = 2500

	# Initialize empty list of k length
	cV_list = np.zeros(k)

	# Iterate for k times in order to retrieve the proper bootstrap sampling
	for i in range(k):

		# Initialize an empty array of length n in order to fill it with 
		# samples for bootstrap resampling
		testArr = np.zeros(n)
		
		# Now select n random numbers from the array. Do this using
		# the numpy random function generator
		for j in range(n):
			
			# Choose a random index coordinate for the sampling
			pos = np.random.randint(n)
			testArr[j] = arr[pos]
		
		# Calculate now the heat capacity from the data on internal energy
		# Do this by using the heat capacity formula
		mean_E = np.mean(testArr)
		meansq_E = np.mean(testArr**2)
		cV = (meansq_E - mean_E**2) / (no_particles * T**2)
		cV_list[i] = cV
	
	# Now calculate the error for the heat capacity from the given data on the
	# retrieved energies	
	mean_cV = np.mean(cV_list)
	meansq_cV = np.mean(cV_list**2)
	err_cV = np.sqrt(meansq_cV - mean_cV**2)
	return err_cV
		
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
					
				# Function for updating the lists for energy and magnetization of the
				# lattice as function of temperature. We start measuring the energy and
				# magnetization when we overreach 100 measurements to lose the initial
				# condition. Moreover, we measure once per 10 measurements, variable 
				# declared at the beginning of the code (lines 25 and26)
				if((dummy_variable % (no_particles * no_sweeps_between) == 0) and (dummy_variable > (no_particles * no_initial_sweeps))):
					time.append(dummy_variable)
					M = LatticeMagnetization(spin_lattice)
					E= LatticeEnergy(spin_lattice)
					magnetizationLattice.append(np.abs(M))
					magnetizationSquared.append((np.abs(M))**2)
					energy.append(E)
					energySquared.append(E**2)
				
				# Exit the function. This happens when we overreach 1000 measurements per
				# temperature. When this happens, perform the Bootstrap sampling to get for
				# the array of energies and magnetization the errors on heat capacity and
				# magnetic susceptibility
				if(dummy_variable > 25250000):
				
					# Perform the bootstrap test
					err_cV = Bootstrap_Sampling(energy)
					err_chi = Bootstrap_Sampling(magnetizationLattice)
					return time, np.mean(magnetizationLattice), np.mean(magnetizationSquared), np.mean(energy), 						np.mean(energySquared), err_cV, err_chi
					
					# Condition for exiting the while loop		
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

def KawasakiUpdate(T):
	dummy_variable = 0
	counter = 0
	
	# Calculate the initial energy of the Ising Model lattice
	E0 = LatticeEnergy(spin_lattice)
	
	# Initialize lists for the measurement and the overall magnetization. These two numpy
	# arrays will be plotted in the end on Matplotlib
	energy = []
	energySquared = []
	time = []
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
				# if(dummy_variable % 1000 == 0):
				# 	image = plt.imshow(spin_lattice, animated=True)
				# 	plt.draw()
				# 	plt.pause(0.01)
				
				# Function for plotting the magnetization of the lattice
				# as function of temperature. The same steps are performed 
				# as for Glauber Dynamics
				if((dummy_variable % (no_particles * no_sweeps_between) == 0) and (dummy_variable > (no_particles * 				no_initial_sweeps))):
					time.append(dummy_variable)
					E= LatticeEnergy(spin_lattice)
					energy.append(E)
					energySquared.append(E**2)
					
				# Exit the function
				if(dummy_variable > 25250000):
				
					# Perform the bootstrap test
					err_cV = Bootstrap_Sampling(energy)
					return time, np.mean(energy), np.mean(energySquared), err_cV		
					counter = 1

# This function represents the iteration through all the temperatures by using the 
# Glauber Dynamics. The function invokes the Glauber Update function, which is only
# performed for one specific temperature
def GlauberDynamics():

	# Initialize the empty arrays for storing values during the measurements.
	# We store values for energy and magnetization lattice, as well as heat
	# capacity and magnetic susceptibility
	temperature = []
	magnetization = []
	susceptibility = []
	energy = []
	heatCapacity = []
	err_heatCapacity = []
	err_susceptibility = []
	count = 0
	
	# Initialize the temperature as 1 and set the while instruction for iterating
	# through the temperatures from 1 to 3 with step of 0.1
	T = 1
	while(T < 3):
		t1 = time.time()
		
		# Return the values when invoking the GlauberUpdate function.
		t, m, msq, E, Esq, err_cV, err_chi = GlauberUpdate(T)
		t2 = time.time()
		
		# Mean magnetic susceptibility and heat capacity have been easy to
		# retrieve directly from the calculations below.
		chi_M = (msq - m**2) / (no_particles * T)
		C_v = (Esq - E**2) / (no_particles * T**2)
		
		# Append the proper values to the lists
		temperature.append(T)
		magnetization.append(m)
		energy.append(E)
		susceptibility.append(chi_M)
		heatCapacity.append(C_v)
		err_heatCapacity.append(err_cV)
		err_susceptibility.append(err_chi)
		
		# Test printing statements
		print("step " + str(count) + " finished with the taken time: " + str((t2-t1)/60) + " minutes")
		print("mean magnetization: " + str(m))
		print("mean lattice energy: " + str(E))
		print("magnetic susceptibility: " + str(chi_M) + "\n")
		print("error on susceptibility: " + str(err_chi))
		print("heat capacity: " + str(C_v))
		print("error on heat capacity: " + str(err_cV) + "\n")
		count += 1
		T += 0.1
	
	# Save all the arrays to a pandas data frame which is afterwards
	# converted to a .csv file
	df = pd.DataFrame({"Temperature" : temperature, "Magnetization" : magnetization, "Lattice Energy" : energy, "Heat Capacity" : 		heatCapacity, "Magnetic Susceptibility" : susceptibility, "Error on heat capacity" : err_heatCapacity, "Error on susceptibility: 		" : err_susceptibility})
	df.to_csv("MCIsingModel.csv")

# Same procedure is applied for the Kawasaki Dynamics	
def KawasakiDynamics():

	# The initial lists for the energy, heat capacity and the errors on
	# the heat capacity are implemented here
	temperature = []
	energy = []
	heatCapacity = []
	err_heatCapacity = []
	
	# Once again, iterathe through the while-instruction from temperature of
	# 1 to the temperature of 3 with step 0.1
	count = 0
	T = 1
	while(T < 3):
		t1 = time.time()
		t, E, Esq, err_cV = KawasakiUpdate(T)
		t2 = time.time()
		C_v = (Esq - E**2) / (no_particles**2 * T**2)
		temperature.append(T)
		energy.append(E)
		heatCapacity.append(C_v)
		err_heatCapacity.append(err_cV)
		print("step " + str(count) + " finished with the taken time: " + str((t2-t1)/60) + " minutes")
		print("mean lattice energy: " + str(E))
		print("heat capacity: " + str(C_v))
		print("error on heat capacity: " + str(err_cV) + "\n")
		count += 1
		T += 0.1

	throw
	df = pd.DataFrame({"Temperature" : temperature, "Lattice Energy" : energy, "Heat Capacity" : heatCapacity, "Error on heat 					capacity" : err_heatCapacity})
	df.to_csv("MCIsingModelKawasaki.csv")

# Now let the user decide whether the Glauber or Kawasaki Dynamics 
# is the prefered choice	
choiceNumber = float(input("Type 0 (Glauber Dynamics) or 1 (Kawasaki Dynamics):"))
if(choiceNumber == 0):
	print("You chose the Glauber Dynamics model!")
	GlauberDynamics()
elif(choiceNumber == 1):
	print("You chose the Kawasaki Dynamics model!")
	KawasakiDynamics()
else:
	raise ValueError("You must choose either 0 or 1!")
