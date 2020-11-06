###################################################
# Plot G(r) using G_r_Ar.txt
############################################

import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore

# Load g(r) data from simulation
# Column 1: distance in meters
# Column 2: number of points in that bin
data = np.loadtxt('G_r_all_distances_O2.txt')

NUM_TIME_STEPS = (200000)/10         # Number of time steps
SIGMA = 0.3006000E-9             

data[:, 0] = data[:, 0]       # Convert distance to units of sigma 
LENGTH = (2.8046445E-9)          # Calculate volume in units of sigma
Volume = LENGTH*LENGTH*LENGTH
NATOMS = 500.0*2.0

bins_list = data[:, 0]
results = data[:, 1]
nbins = results.shape[0]
g_of_r = np.zeros(nbins)
scale_factor = np.zeros(500)

# Calculate scale factor
step = bins_list[1]-bins_list[0]
step = 0.05E-10
for i in range(1, nbins):
    expected = NATOMS*np.pi*4.0/Volume*((i*step+step)**3-(i*step)**3 )/3
    scale_factor[i] = expected*NUM_TIME_STEPS*NATOMS
    #scale_factor[i] = NATOMS*NATOMS*4.0*np.pi*step*(i*step)**2
   
print(step)
# Scale g(r)
for i in range(1, nbins):
   # g_of_r[i] = (results[i])*Volume/(scale_factor[i])/NUM_TIME_STEPS
    g_of_r[i] = (results[i])/scale_factor[i]

    bins_list[i] = step*i/1.0E-9


# Plot
plt.rcParams.update({'font.size': 20})
plt.figure(figsize=(15,8))
plt.plot(bins_list[0:nbins], g_of_r[0:nbins])
plt.plot([0, bins_list[nbins-1]], [1, 1], 'k')
plt.xlabel('$r$ (nm)')
plt.ylabel('$g(r)$')
plt.axis([0, 1.25, 0, 3])
plt.minorticks_on()
plt.savefig('O2_Gr.png', bbox_inches='tight')
g1 = np.concatenate((bins_list, g_of_r))
np.savetxt('O2_gr_data.txt', g1)
plt.show()
