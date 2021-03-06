'''

'''
from trojan import *
# import matplotlib.pyplot as plt
# import numpy as np
# import math
# from trojan import (
# 	Two_Body_System,
# 	measure_period,
# 	t_span, 
# 	N, 
# 	deltas, 
# 	rand_deltas,
# 	ODE,
# 	Point,
# 	Conditions
# 	)
system = Two_Body_System(m_j, m_s)

# rand_vels = np.loadtxt('data/rand_vels.txt')
# rand_deltas = np.loadtxt('data/rand_deltas.txt')
# deltas = get_deltas(asteroids)
# delta_vs = get_delta_vs(asteroids)
dxs = get_ds(min_delta, max_delta, asteroids)
dys = get_ds(min_delta, max_delta, asteroids)

# def make_plots_ivp():

# 	plt.plot(ys[0][0], ys[0][2], label = "Jupiter")
# 	plt.plot(ys[0][4], ys[0][6], label = "Sun")
# 	# plt.plot(y[8], y[10], label = "Asteroid")
# 	plt.legend()
# 	plt.show()


unpert_ode = ODE(np.loadtxt('data/unperturbed.txt'))

# solved_deltas = []
# for i, delta in enumerate(deltas):
# 	for j, delta_v in enumerate(delta_vs):
# 		solved_delta = ODE(np.loadtxt(f'data/delta_{i}_{j}.txt'))
# 		solved_deltas.append(solved_delta)

# solved_rand_deltas = []
# for i, rand_delta in enumerate(rand_deltas):
# 	for j, rand_vel in enumerate(rand_vels):
# 		solved_rand_delta = ODE(np.loadtxt(f'data/rand_delta_{i}_{j}.txt'))
# 		solved_rand_deltas.append(solved_rand_delta)

single_deltas = []
for i, dy in enumerate(dys):
	single_delta = ODE(np.loadtxt(f'data/single_delta{i}.txt'))
	single_deltas.append(single_delta)

# solved_grids = []
# for i, dx in enumerate(dxs):
# 	for j, dy in enumerate(dys):
# 		solved_grid = ODE(np.loadtxt(f'data/position_grid_{i}_{j}.txt'))
# 		solved_grids.append(solved_grid)





# deviation(unperturbed, perturbed)
system.planar_plot(-8, 8, -8, 8, single_deltas, unpert_ode)
# stationary_initial_position(unpert_ode, solved_deltas)
# stationary_initial_position(unpert_ode, solved_rand_deltas)
corotated_deviation(unpert_ode, single_deltas, "Motion in Co-Rotating Frame")
# corotated_deviation(unpert_ode, solved_rand_deltas, "Motion in Co-Rotating Frame (random displacements)")
# make_plots_ode(unpert_ode, single_deltas)

plt.show()


