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

rand_vels = np.loadtxt('data/rand_vels.txt')
rand_deltas = np.loadtxt('data/rand_deltas.txt')
deltas = get_deltas(asteroids)
delta_vs = get_delta_vs(asteroids)
dxs = get_ds(asteroids)
dys = get_ds(asteroids)
dv_xs = get_ds(asteroids)
dv_ys = get_ds(asteroids)

unpert_ode = ODE(np.loadtxt('data/unperturbed.txt'))

solved_rand_deltas = []
for i, rand_delta in enumerate(rand_deltas):
	for j, rand_vel in enumerate(rand_vels):
		solved_rand_delta = ODE(np.loadtxt(f'data/rand_delta_{i}_{j}.txt'))
		solved_rand_deltas.append(solved_rand_delta)

solved_deltas = []
for i, delta in enumerate(deltas):
	for j, delta_v in enumerate(delta_vs):
		solved_delta = ODE(np.loadtxt(f'data/delta_{i}_{j}.txt'))
		solved_deltas.append(solved_delta)

solved_grids = []
for i, dx in enumerate(dxs):
	for j, dy in enumerate(dys):
		solved_grid = ODE(np.loadtxt(f'data/grid_{i}_{j}.txt'))
		solved_grids.append(solved_grid)

solved_vs = []
for i, dv_x in enumerate(dv_xs):
	for j, dv_y in enumerate(dv_ys):
		solved_v = ODE(np.loadtxt(f'data/v_{i}_{j}.txt'))
		solved_vs.append(solved_v)



# def maximum_deviation():
# 	max_devs = []
# 	r_primed_unpert = perform_rotation_ode(unpert_ode)

# 	for i, solved_rand_delta in enumerate(solved_rand_deltas):
# 		print(f"solved_rand_delta: {solved_rand_delta}")
# 		r_primed_rand_delta = perform_rotation_ode(solved_rand_delta)
# 		vector_devs = np.array(r_primed_rand_delta) - np.array(r_primed_unpert)
# 		print(f"vector devs: {vector_devs}")
# 		devs = np.linalg.norm(vector_devs, axis = 1)
# 		max_dev = np.amax(devs)
# 		print(f"the maximum deviation for {i} is {max_dev}")

# 	max_devs.append(max_dev)
# 	print(f"maxdevs: {max_devs}")
# 	plt.plot(rand_deltas, max_devs)
# 	plt.title("Maximum deviation from L4 versus initial displacement")

# 	# return(max_devs)



# maximum_deviation(solved_deltas, deltas)
# heatmap(solved_deltas, unpert_ode, "Change in velocity in x direction /arb units", "Displacement from L4 in y direction / arb units")
heatmap(solved_grids, unpert_ode, "dx /arb units", "dy / arb units")
heatmap(solved_vs, unpert_ode, "dv_x /arb units", "dv_y / arb units")


# plt.plot(mod_deltas, max_mod_delta_rs)
# plt.title("Max deviation vs initial displacement from L4")

plt.show()



