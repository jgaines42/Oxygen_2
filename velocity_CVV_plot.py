###################################################
# Plot CVV autocorrelation
############################################
import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import scipy as scipy
from scipy import optimize  # type: ignore
import math


# Define to fit exponential
def exponential(x, a, k):
    return a * np.exp(x * k)


# Load file
# Velocity_Ar.txt is the autocorrelation data from the run
data = np.loadtxt('Velocity_O2.txt')

average_auto = data
DT = 4.0E-15 / (1E-12)       # Time step in ps
DT_fs = 4.0E-15            # Time step in s
NUMBER_ATOMS = 500          # Number of atoms
NUMBER_TIME = 300000         # Number of time steps
EQUILIBRIUM_START = 0       # When equlibirum starts (in the NUMBER_TIME)
number_steps = 1000
NUMBER_STARTS = NUMBER_TIME - EQUILIBRIUM_START - number_steps
print(NUMBER_STARTS)

# Normalilze data by the number of time steps in it and the number of atoms
average_auto = np.divide(average_auto, NUMBER_STARTS)
average_auto = np.divide(average_auto, NUMBER_ATOMS)

# Make time data
all_steps = np.zeros([len(data)])
for i in range(0, len(average_auto)):
    all_steps[i] = i * DT         # in ps

# Plot CVV(t)
plt.figure(figsize=(15, 8))
plt.rcParams.update({'font.size': 20})
plt.xlabel('time (ps)')
plt.ylabel('$C_{vv}(t), (\mathrm{m^2/s^2})$')
plt.plot(all_steps, average_auto, 'k')
plt.plot([0, all_steps[number_steps - 1]], [0, 0])
plt.axis([0, 4, -10000, 61000])
plt.savefig('Ar_CVV.png', bbox_inches='tight')
plt.show()

print(average_auto[0])

# Plot normalized by dt = 0
plt.figure(figsize=(15, 8))
plt.plot(all_steps, average_auto / average_auto[0], 'k')
plt.xlabel('time (ps)')
plt.ylabel('$C_{vv}(t), (\mathrm{m^2/s^2})$')
plt.axis([0, 1.5, -0.2, 4])
plt.show()
plt.savefig('O2_CVV_normalized.png', bbox_inches='tight')

np.savetxt('O2_CVV_data.txt', average_auto)


# Use time in seconds to fit data
all_steps = np.zeros([len(data)])
for i in range(0, len(average_auto)):
    all_steps[i] = i * DT_fs         # in s

# Fit tail of data to exponential
print(all_steps[250])
popt_exponential, pcov_exponential = scipy.optimize.curve_fit(exponential, all_steps[250:number_steps],
                                                              average_auto[250:number_steps], p0=[-20000, -3E12])
perr_exponential = np.sqrt(np.diag(pcov_exponential))
print("pre-exponential factor = %0.2f (+/-) %0.2f" % (popt_exponential[0], perr_exponential[0]))
print("rate constant = %0.2f (+/-) %0.2f" % (popt_exponential[1], perr_exponential[1]))

# Plot exponential equation compared to data
plt.figure(figsize=(15, 8))
plt.plot(all_steps[100:number_steps], exponential(all_steps[100:number_steps],
         popt_exponential[0], popt_exponential[1]), 'k--')
plt.plot(all_steps[100:number_steps], average_auto[100:number_steps])
plt.savefig('cvv_fit.png', bbox_inches='tight')

plt.show()

# Extract k and a and convert k to seconds
k = popt_exponential[1]  # 1/s
a = popt_exponential[0]  # m^2


# Integrate from all_steps[150] to infinity
# intergral -> 1/k *exp(kx)

end_step = all_steps[number_steps - 1] + DT_fs
integral_15 = np.sum(average_auto[0:number_steps - 1]) * DT_fs     # m^2/s
integral_15p = 0.0 - (1.0 / k * a * math.exp(k * end_step))        # m^2/ps

print(integral_15 * (100 * 100))
print(integral_15p * (100 * 100))

# Calculate diffusion coefficient
D = (integral_15 + integral_15p) / 3.0  # m^2/s
D = (D) * (100 * 100)                   # cm^2/s
print(D)
