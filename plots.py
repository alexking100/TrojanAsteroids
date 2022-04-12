'''

'''

import matplotlib.pyplot as plt
import numpy as np
import math
from trojan import Two_Body_System
from trojan import measure_period
from trojan import t_span
from trojan import N
from trojan import deltas, rand_deltas

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



def make_plots_ivp():

	plt.plot(ys[0][0], ys[0][2], label = "Jupiter")
	plt.plot(ys[0][4], ys[0][6], label = "Sun")
	# plt.plot(y[8], y[10], label = "Asteroid")
	plt.legend()
	plt.show()


def make_plots_ode():

	plt.plot(unpert_ode.r_j.x, unpert_ode.r_j.y, label = "Jupiter")
	plt.plot(unpert_ode.r_s.x, unpert_ode.r_s.y, label = "Sun")
	plt.plot(unpert_ode.r_a.x, unpert_ode.r_a.y, label = "Unperturbed Asteroid")

	plt.plot(solved_deltas[0].r_a.x, solved_deltas[0].r_a.y, label = "Perturbed Asteroid")
	plt.legend()
	plt.show()


def angle(t_eval):

	T = measure_period(unpert_ode)

	omega = 2*math.pi / T
	thetas = np.zeros(len(t_eval))

	for i, t in enumerate(t_eval):
		thetas[i] = omega*t

	return thetas



def perform_rotation_ode(ode):
	'''
	rotate the x-y frame at the same rate as the eqm asteroid
	returns two vectors in a tuple, each length of t_eval and containing x/y coordinates in the co-rotated frame
	'''
	t_eval = np.linspace(t_span[0], t_span[1], N)
	thetas = angle(t_eval)

# for asteroid
	x_primed = np.cos(thetas) * ode.r_a.x - np.sin(thetas) * ode.r_a.y
	y_primed = np.sin(thetas) * ode.r_a.x + np.cos(thetas) * ode.r_a.y

# for jupiter
	# x_primed = np.cos(thetas) * ode.r_j.x - np.sin(thetas) * ode.r_j.y
	# y_primed = np.sin(thetas) * ode.r_j.x + np.cos(thetas) * ode.r_j.y	

# for sun
	# x_primed = np.cos(thetas) * ode.r_s.x - np.sin(thetas) * ode.r_s.y
	# y_primed = np.sin(thetas) * ode.r_s.x + np.cos(thetas) * ode.r_s.y	

	r_primed = np.array([x_primed, y_primed])

	return r_primed


def corotated_deviation(unpert_ode, solved_deltas):

	print(len(solved_deltas))

	# diff_xs = np.zeros(len(delta_odes))
	# diff_ys = np.zeros(len(delta_odes))

	# r_primed_eqm, r_primed_pert = perform_rotation(ys)
	r_primed_unpert = perform_rotation_ode(unpert_ode)
	plt.plot(r_primed_unpert[0], r_primed_unpert[1], label = "L4")

	for i, solved_delta in enumerate(solved_deltas):

		r_primed_delta = perform_rotation_ode(solved_delta)

		# diff_xs[i] = r_primed_delta[0] - r_primed_unpert[0]
		# diff_ys[i] = r_primed_delta[1] - r_primed_unpert[1]

		plt.plot(r_primed_delta[0], r_primed_delta[1], label = f"delta = {float(round(np.linalg.norm(deltas[i]), 3))}")

	plt.legend()
	plt.title("Co-Rotated Motion of Asteroids")
	plt.show()

	# diff_ast_x = r_primed_unpert[0] - unpert_ode.r_a.x[0]
	# diff_ast_y = r_primed_unpert[1] - unpert_ode.r_a.y[0]
	# print(unpert_ode.r_a.x[0], unpert_ode.r_a.y[0])

	# plt.plot(unpert_ode.r_a.x, unpert_ode.r_a.y, label = "Asteroid in rest frame")
	# plt.plot(r_primed_unpert[0], r_primed_unpert[1], label = "Asteroid in rotated frame")
	# plt.legend()

	# fig, ax = plt.subplots()
	# plt.plot(diff_x, diff_y, label = 'Deviation of asteroid') # difference between the asteroid's position and the initial lagrange point



    


# y = np.loadtxt('orbitdata.txt')

unperturbed = np.loadtxt('unperturbed.txt')
unpert_ode = ODE(unperturbed)

solved_deltas = []
for i, delta in enumerate(deltas):
	solved_delta = ODE(np.loadtxt(f'delta_{i}.txt'))
	solved_deltas.append(solved_delta)

solved_rand_deltas = []
for i, rand_delta in enumerate(rand_deltas):
	solved_rand_delta = ODE(np.loadtxt(f'rand_delta_{i}.txt'))
	solved_rand_deltas.append(solved_rand_delta)




# deviation(unperturbed, perturbed)


# corotated_deviation(unpert_ode, solved_deltas)
# corotated_deviation(unpert_ode, solved_rand_deltas)
# make_plots_ode()


