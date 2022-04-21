'''
This script will test the stability of orbits with a 
random initial displacement from L4 and a random velocity (deviated from equilibrium velocity)

'''
from trojan import *
# import matplotlib.pyplot as plt
# import numpy as np
# import math
# from trojan import deltas, rand_deltas, asteroids
# from plots import perform_rotation_ode
# from two_body_system import Two_Body_System, Point, ODE
system = Two_Body_System(m_j, m_s)

# rand_vels = np.loadtxt('data/rand_vels.txt')
# rand_deltas = np.loadtxt('data/rand_deltas.txt')
deltas = get_deltas(asteroids)
delta_vs = get_delta_vs(asteroids)
dxs = get_ds(-max_delta, max_delta, asteroids)
dys = get_ds(-max_delta, max_delta, asteroids)
dv_xs = get_ds(-max_delta_v, max_delta_v, asteroids)
dv_ys = get_ds(-max_delta_v, max_delta_v, asteroids)

unpert_ode = ODE(np.loadtxt('data/unperturbed.txt'))

# solved_rand_deltas = []
# for i, rand_delta in enumerate(rand_deltas):
# 	for j, rand_vel in enumerate(rand_vels):
# 		solved_rand_delta = ODE(np.loadtxt(f'data/rand_delta_{i}_{j}.txt'))
# 		solved_rand_deltas.append(solved_rand_delta)

# solved_deltas = []
# for i, delta in enumerate(deltas):
# 	for j, delta_v in enumerate(delta_vs):
# 		solved_delta = ODE(np.loadtxt(f'data/delta_{i}_{j}.txt'))
# 		solved_deltas.append(solved_delta)

# solved_grids = []
# for i, dx in enumerate(dxs):
# 	for j, dy in enumerate(dys):
# 		solved_grid = ODE(np.loadtxt(f'data/position_grid_{j}_{i}.txt'))
# 		solved_grids.append(solved_grid)

# solved_vs = []
# for i, dv_x in enumerate(dv_xs):
# 	for j, dv_y in enumerate(dv_ys):
# 		solved_v = ODE(np.loadtxt(f'data/velocity_grid_{i}_{j}.txt'))
# 		solved_vs.append(solved_v)

solved_dys = []
for i, dy in enumerate(dys):
	solved_dy = ODE(np.loadtxt(f'data/solved_dys_{i}.txt'))
	solved_dys.append(solved_dy)

maximum_deviation(solved_dys, dys, unpert_ode, direction = "y")
maximum_deviation(solved_dys, dys, unpert_ode, direction = "y", cap = True)

# def maximum_deviation(solved_deltas):
# 	max_devs = []
# 	r_primed_unpert = perform_rotation_ode(unpert_ode, unpert_ode)

# 	for i, solved_delta in enumerate(solved_deltas):
# 		r_primed_delta = perform_rotation_ode(unpert_ode, solved_delta)
# 		vector_devs = np.array(r_primed_delta) - np.array(r_primed_unpert)
# 		devs = np.linalg.norm(vector_devs, axis = 1)
# 		max_dev = np.amax(devs)
# 		max_devs.append(max_dev)

# 	fig, ax = plt.subplots()
# 	ax.plot(np.linalg.norm(deltas), max_devs)
# 	ax.set_ylim(0, 5)
# 	ax.set_title("Maximum deviation from L4 versus initial displacement")

# 	# return(max_devs)


# maximum_deviation(solved_deltas, deltas)
# heatmap(solved_deltas, unpert_ode, "Change in velocity in x direction /arb units", "Displacement from L4 in y direction / arb units")
# stationary_initial_position(unpert_ode, solved_deltas, "x /AU", "y /AU", "Initial positions of Asteroids")
# stationary_initial_position(unpert_ode, solved_vs, "v_x", "v_y", "Initial velocities of Asteroids")

# maximum_deviation(solved_deltas)

# system.planar_plot(-8, 8, -8, 8, solved_grids, unpert_ode) # plot for initial displacements of asteroids for the grid
# system.planar_plot(3, 7, 1, 5, solved_grids, unpert_ode) # plot for initial displacements of asteroids for the grid
# system.planar_plot(unpert_ode.r_a.x[0] - max_delta, unpert_ode.r_a.x[0] + max_delta, unpert_ode.r_a.y[0] - max_delta, unpert_ode.r_a.y[0] + max_delta, solved_grids, unpert_ode)

# heatmap(solved_grids, unpert_ode, "dx /arb units", "dy / arb units")
# heatmap(solved_vs, unpert_ode, "dv_x /arb units", "dv_y / arb units")


# plt.plot(mod_deltas, max_mod_delta_rs)
# plt.title("Max deviation vs initial displacement from L4")

plt.show()



