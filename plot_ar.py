###################################################
# Plot energy, temperature and speed distributions
############################################

import math  # type: ignore
import numpy as np  # type: ignore
import matplotlib.pyplot as plt   # type: ignore

##############################################
# Plot energy
##############################################

NVE_START = 0
NUM_RUNS = 100000
MASS = 0.016 * 2.0 / 6.02E23              # Atomic mass
Timestep = 4.0E-15
data = np.loadtxt('Energies_O2.txt')
data[:, 0] = np.multiply(data[:, 0], 1.0E12)    # Change time to ps

# Format plot
plt.figure(figsize=(15, 8))
plt.rcParams.update({'font.size': 20})

# Plot all energy terms in fJ
plt.plot(data[:, 0], np.multiply(data[:, 1], 1.0E15), label='Potential')
plt.plot(data[:, 0], np.multiply(data[:, 2], 1.0E15), label='Kinetic')
plt.plot(data[:, 0], np.multiply(data[:, 3], 1.0E15), label='Total')
plt.xlabel('time (ps)')
plt.ylabel('Energy (fJ)')
plt.legend(loc='center right')
plt.savefig('O2_Energy.png', bbox_inches='tight')
plt.show()


##############################################
# Plot Temperature
##############################################

dt_NVE = (Timestep) * 1.0E12          # Timestep in NVE in ps
data3 = np.loadtxt('Temp.txt')   # NVE Temp

data4 = np.arange(data3.shape[0])
data4 = data4 * dt_NVE                # Time in NVE

# Plot NVE temperature
plt.figure(figsize=(15, 8))
plt.plot(data4, data3, 'k')
plt.xlabel('time (ps)')
plt.ylabel('Temperature (K)')
plt.plot([0, data4[data4.shape[0] - 1]], [np.mean(data3[NVE_START:NUM_RUNS]), np.mean(data3[NVE_START:NUM_RUNS])], 'r')
plt.plot([0, data4[data4.shape[0] - 1]], [77, 77], '--r')

plt.savefig('O2_Temp_NVE.png', bbox_inches='tight')
plt.show()

# Print average and std of temp of NVE
print(np.mean(data3[NVE_START:NUM_RUNS]))
print(np.std(data3[NVE_START:NUM_RUNS]))

# Calculate average temp in NVE
mean_temp = np.mean(data3[NVE_START:NUM_RUNS])
massO = 0.016 * 2.0 / (6.02E23)
eq_v = (3.0 * (1.38064852E-23) * mean_temp / massO)
print(eq_v)

##############################################
# Plot Speed distribution
##############################################

data = np.loadtxt('Vel_distribution_O2.txt')

T = mean_temp

KB = 1.38064852E-23
EPSILON = 120.0 * KB
SCALE_FACTOR = math.sqrt(math.pow(MASS / (2.0 * np.pi * KB * T), 3))

mb = np.zeros([800, 2])                     # maxwell boltzman speed distribution
TARGET_V = math.sqrt(3.0 * KB * T / MASS)         # Initial velocity of simulation
print(TARGET_V)

# sample maxwell boltzman from 0 to 600 m/s
STEP = 1
for i in range(1, 800):
    this_v = STEP * i                         # Current velocity
    mb[i, 0] = (math.exp(-MASS * (this_v) * (this_v) / (2 * KB * T)) * SCALE_FACTOR * 4 * np.pi * this_v * this_v)
    mb[i, 0] = mb[i, 0]                     # P(v) at this time point
    mb[i, 1] = this_v

# Plot data
plt.figure(figsize=(15, 8))
plt.rcParams.update({'font.size': 20})
plt.plot(mb[:, 1], mb[:, 0], 'k')
plt.plot(mb[:, 1], data / np.sum(data))
plt.xlabel('speed (m/s)')
plt.ylabel('P(speed)  (s/m)')
plt.axis([0, 800, 0, 0.005])
plt.savefig('O2_vel_dist.png', bbox_inches='tight')
plt.show()
