import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
data = np.loadtxt('Velocity_Vector.txt')

auto_coor = np.zeros([864, 190])

average_auto = np.zeros([50])
DT = 10.0E-15/(1E-12)
NUMBER_STARTS = 100
EQUILIBRIUM_START = 2500
all_steps = np.zeros(50)

for step_loop in range(0, 50):
    all_steps[step_loop] = step_loop*10.0*DT
    step_size = step_loop*10.0

    counter = 0
    for time_start in range(0, NUMBER_STARTS):

        # velocities sare saved every 10 time steps
        counter = counter + 1
        t1 = time_start*10+EQUILIBRIUM_START
        t2 = time_start*10+EQUILIBRIUM_START+step_size

        ind0 = data[:, 0] == t1
        ind1 = data[:, 0] == t2
        this_1 = data[ind0, 1:4]
        this_2 = data[ind1, 1:4]
        mult1 = np.multiply(this_1, this_2)
        average_auto[step_loop] = average_auto[step_loop]+np.sum(mult1)


average_velocity = np.multiply(data[:, 1], data[:, 1])
average_velocity = average_velocity + np.multiply(data[:, 2], data[:, 2])
average_velocity = average_velocity + np.multiply(data[:, 3], data[:, 3])
average_velocity = np.sum(average_velocity)/data.shape[0]

average_auto = np.divide(average_auto, NUMBER_STARTS)
average_auto = np.divide(average_auto, 864.0)
average_auto = np.divide(average_auto, average_velocity)

print(average_auto[0])
plt.rcParams.update({'font.size': 20})
plt.xlabel('time (ps)')
plt.ylabel('$C_{vv}(t)$')
plt.plot(all_steps[0:49], average_auto[0:49], 'k')
plt.plot([0, all_steps[49]], [0, 0])
plt.axis([0, 5, -.1, 1])

plt.show()
