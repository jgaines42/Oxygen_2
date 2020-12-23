#############################################
# Plot mean squared displacement using data generated in fortran code
#############################################
import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import scipy as scipy
from scipy import optimize


def linear(x, a, b):
    return a * x + b


data = np.loadtxt('Mean_squared_displacement.txt')

DT = 4.0E-15
natoms = 500
num_runs = 300000
num_msd = 3000
all_steps = np.zeros(num_msd)

# Normalize data
for step_loop in range(0, num_msd):
    all_steps[step_loop] = step_loop * DT * 1E12                 				# time in ps
    data[step_loop] = data[step_loop] / natoms / (num_runs - num_msd) * (1E20)  # msd in Angstroms^2

fit_start = 1500
print(all_steps[fit_start])
print(all_steps[2999])

# Fit linear section to calculate diffusion coeficient
popt_linear, pcov_linear = scipy.optimize.curve_fit(linear, all_steps[fit_start:3000], data[fit_start:3000],
                                                    p0=[8000.0, 0])
perr_linear = np.sqrt(np.diag(pcov_linear))
print("slope = %0.2f (+/-) %0.2f" % (popt_linear[0], perr_linear[0]))
print("intercept= %0.2f (+/-) %0.2f" % (popt_linear[1], perr_linear[1]))

m = popt_linear[0]
b = popt_linear[1]
diffusion = popt_linear[0] / 6.0  # diffusion in A^2/ps
print(diffusion * 100.0 * 100.0 * (1E12) / (1E20))  # Difusion in cm^2/s


plt.figure(figsize=(15, 8))
plt.rcParams.update({'font.size': 20})
plt.plot(all_steps, data, 'k', label='Data')
plt.plot([0, all_steps[num_msd - 1]], [b, b + all_steps[num_msd - 1] * m], '--b', label='Fit')
plt.plot([all_steps[fit_start], all_steps[2999]], [data[fit_start], data[2999]], 'x')
plt.xlabel('Time (ps)')
plt.ylabel('Mean squared displacement ($\AA ^2$)')
plt.legend(loc='center right')
plt.minorticks_on()
plt.savefig('O2_msd.png', bbox_inches='tight')

plt.show()
