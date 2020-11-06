import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
data = np.loadtxt('Velocity_Vector.txt')
print(data[0, :])

DT = 10.0E-15
NUMBER_STARTS = 100
EQUILIBRIUM_START = 2500
all_steps = np.zeros(30)
all_msd = np.zeros([NUMBER_STARTS, 30])
for time_start in range(0, NUMBER_STARTS):
    all_coord = np.zeros([864, 3])
    for step_loop in range(0, 30):
        all_steps[step_loop] = step_loop*10.0*DT*10**12
        step_size = step_loop*10.0

        t2 = time_start*10+EQUILIBRIUM_START+step_size
        ind1 = data[:, 0] == t2
        this_2 = data[ind1, 1:4]

        all_coord = all_coord + this_2*DT*10.0
        ms1 = (all_coord[:, 0]*all_coord[:, 0]
               + all_coord[:, 1]*all_coord[:, 1]
               + all_coord[:, 2]*all_coord[:, 2])
        all_msd[time_start, step_loop] = np.sum(ms1)/(864.0)*10**20

plt.plot(all_steps, all_msd[0, :])
plt.plot(all_steps, all_msd[NUMBER_STARTS-1, :])

x = np.mean(all_msd, axis=0)
plt.plot(all_steps, x, 'k')

diffusion = (x[29]-x[19])/(all_steps[29]-all_steps[19])/6.0
print(diffusion)

plt.axis([0, 3, 0, 5])
plt.xlabel('Time (ps)')
plt.ylabel('Mean squared displacement ($\AA ^2$)')
plt.show()
