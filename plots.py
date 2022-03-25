import matplotlib.pyplot as plt
import numpy as np

y = np.loadtxt('orbitdata.txt')


# print(y)
# print(y[9], y[11])
# y = [x_j, vx_j, y_j, vy_j, x_s, vx_s, y_s, vy_s, x_a, vx_a, y_a, vy_a]

plt.plot(y[0], y[2], label = "Jupiter")
plt.plot(y[4], y[6], label = "Sun")
plt.plot(y[8], y[10], label = "Asteroid")
plt.legend()
plt.show()
