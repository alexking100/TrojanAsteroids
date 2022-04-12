'''
This script will test the stability of orbits with a 
random initial displacement from L4 and a random velocity

'''

import matplotlib.pyplot as plt
import numpy as np
import math
# from trojan import Two_Body_System
# from trojan import measure_period
# from trojan import t_span
# from trojan import N
from trojan import deltas, rand_deltas, asteroids
from plots import perform_rotation_ode


class Point:
	def __init__(self, r):
		self.x = r[0]
		self.y = r[1]

class ODE:
	def __init__(self, ys):
		self.ode = ys

		self.r_j = Point([ys[:, 0], ys[:, 2]])
		self.v_j = Point([ys[:, 1], ys[:, 3]])

		self.r_s = Point([ys[:, 4], ys[:, 6]])
		self.v_s = Point([ys[:, 5], ys[:, 7]])

		self.r_a = Point([ys[:, 8], ys[:, 10]])
		self.v_a = Point([ys[:, 9], ys[:, 11]])


unperturbed = np.loadtxt('unperturbed.txt')
unpert_ode = ODE(unperturbed)

solved_rand_deltas = []
for i, rand_delta in enumerate(rand_deltas):
	solved_rand_delta = ODE(np.loadtxt(f'rand_delta_{i}.txt'))
	solved_rand_deltas.append(solved_rand_delta)

def maximum_deviation():

	max_mod_delta_rs = []

	for r in solved_rand_deltas:
		r_primed = perform_rotation_ode(r)
		r_primed_unpert = perform_rotation_ode(unpert_ode)
		delta_r = r_primed - r_primed_unpert
		# print(delta_r)
		mod_delta_r = np.linalg.norm(delta_r, axis = 0)
		# print(mod_delta_r)
		# print(f"shape is {np.shape(mod_delta_r)}")
		max_mod_delta_r = np.amax(mod_delta_r)
		max_mod_delta_rs.append(max_mod_delta_r)

	mod_deltas = np.linalg.norm(rand_deltas, axis = 1)
	# print(max_mod_delta_rs)
	# print(mod_deltas)

	plt.plot(mod_deltas, max_mod_delta_rs, marker = 'o')
	plt.title("Max deviation vs initial displacement from L4")
	plt.show()
	# return (mod_deltas, max_mod_delta_rs)


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


maximum_deviation()

# plt.plot(mod_deltas, max_mod_delta_rs)
# plt.title("Max deviation vs initial displacement from L4")
# plt.show()



