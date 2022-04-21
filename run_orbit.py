from trojan import *

# deltas = get_deltas(asteroids)
# delta_vs = get_delta_vs(asteroids)
# rand_deltas = get_rand_deltas(asteroids)
# rand_vels = get_rand_vels(asteroids)
dxs = get_ds(-max_delta, max_delta, asteroids)
dys = get_ds(-max_delta, max_delta, asteroids)
dv_xs = get_ds(-max_delta_v, max_delta_v, asteroids)
dv_ys = get_ds(-max_delta_v, max_delta_v, asteroids)

'''
Cannot run the random functions again in another file as they will be different
Save them to .txt and call the data in another file
'''
# np.savetxt('data/rand_vels.txt', rand_vels)
# np.savetxt('data/rand_deltas.txt', rand_deltas)


system = Two_Body_System(m_j, m_s)

unperturbed = system.interact_ode(t_span, [0,0], [0,0], v_eqm=True, dv=False)
# print(unperturbed)
# print(np.shape(unperturbed))
np.savetxt('data/unperturbed.txt', unperturbed)

def solver_methods():
    unperturbed_ivp_RK45 = system.interact_ivp(t_span, t, 'RK45', [0,0], [0,0], eqm = True).y
    unperturbed_ivp_RK45 = np.transpose(unperturbed_ivp_RK45)

    unperturbed_ivp_RK23 = system.interact_ivp(t_span, t, 'RK23', [0,0], [0,0], eqm = True).y
    unperturbed_ivp_RK23 = np.transpose(unperturbed_ivp_RK23)

    unperturbed_ivp_DOP853 = system.interact_ivp(t_span, t, 'DOP853', [0,0], [0,0], eqm = True).y
    unperturbed_ivp_DOP853 = np.transpose(unperturbed_ivp_DOP853)

    unperturbed_ivp_Radau = system.interact_ivp(t_span, t, 'Radau', [0,0], [0,0], eqm = True).y
    unperturbed_ivp_Radau = np.transpose(unperturbed_ivp_Radau)

    unperturbed_ivp_BDF = system.interact_ivp(t_span, t, 'BDF', [0,0], [0,0], eqm = True).y
    unperturbed_ivp_BDF = np.transpose(unperturbed_ivp_BDF)

    unperturbed_ivp_LSODA = system.interact_ivp(t_span, t, 'LSODA', [0,0], [0,0], eqm = True).y
    unperturbed_ivp_LSODA = np.transpose(unperturbed_ivp_LSODA)

    print(np.shape(unperturbed_ivp_RK45)) # for some reason, this method overrides the t and plots only 70 points not 10,000
    print(np.shape(unperturbed_ivp_RK23))
    print(np.shape(unperturbed_ivp_Radau))
    print(np.shape(unperturbed_ivp_DOP853))
    print(np.shape(unperturbed_ivp_BDF))
    print(np.shape(unperturbed_ivp_LSODA))


    np.savetxt('data/solver/unperturbed_ivp_RK45.txt', unperturbed_ivp_RK45)
    np.savetxt('data/solver/unperturbed_ivp_RK23.txt', unperturbed_ivp_RK23)
    np.savetxt('data/solver/unperturbed_ivp_DOP853.txt', unperturbed_ivp_DOP853)
    np.savetxt('data/solver/unperturbed_ivp_Radau.txt', unperturbed_ivp_Radau)
    np.savetxt('data/solver/unperturbed_ivp_BDF.txt', unperturbed_ivp_BDF)
    np.savetxt('data/solver/unperturbed_ivp_LSODA.txt', unperturbed_ivp_LSODA)

# solver_methods()


def get_velocity(dx, dy):
    '''
    gets velocity for a point NOT at the lagrange point
    '''

    initial_conditions = circ_orbit_conditions_new(R, m_j, m_s, [0,0], [0,0], [0,0], v_eqm = True, dv = False)
    initial_conds = Conditions(initial_conditions)

    # r is the displacement from the origin of the asteroids, which are displaced from L4 by (dx, dy)
    r = np.array([initial_conds.r_a.x, initial_conds.r_a.y]) + np.array([dx, dy])
    omega = 2 * math.pi / T

    mod_v = np.linalg.norm(r)*omega
    theta = np.arctan(r[0] / r[1])

    v = np.array([mod_v*np.cos(theta), - mod_v*np.sin(theta)])
    print(mod_v)
    # print(theta)
    return v

# get data for a line of asteroids in the y direction
for i, dy in enumerate(dys):
    solved_dys = system.interact_ode(t_span, [0,dy], get_velocity(0, dy), v_eqm = False, dv=False)
    np.savetxt(f'data/solved_dys_{i}.txt', solved_dys)







# get data for a grid of displaced asteroids in the x-y plane about L4
# for i, dx in enumerate(dxs):
#     for j, dy in enumerate(dys):
#         solved_grid = system.interact_ode(t_span, [dx,dy], get_velocity(dx, dy), v_eqm = False, dv=False)
#         np.savetxt(f'data/position_grid_{i}_{j}.txt', solved_grid)

# get data for a grid of velocities to add to the equilibrium velocity, all asteroids are positioned at L4
# for i, dv_x in enumerate(dv_xs):
#     for j, dv_y in enumerate(dv_ys):
#         solved_v = system.interact_ode(t_span, [0,0], [dv_x, dv_y], v_eqm=True, dv=True)
#         np.savetxt(f'data/velocity_grid_{i}_{j}.txt', solved_v)

# get data for formulaically displaced asteroids in y, v_x space
# solved_deltas = []
# for i, delta in enumerate(deltas):
#     # print(delta)
#     for j, delta_v in enumerate(delta_vs): 
#         solved_delta = system.interact_ode(t_span, delta, delta_v, v_eqm=True, dv=True)
#         solved_deltas.append(solved_delta)
#         np.savetxt(f'data/delta_{i}_{j}.txt', solved_delta)

# get data for randomly displaced asteroids (each delta value has every random velocity attached)
# solved_rand_deltas is of size (asteroids*asteroids)
# solved_rand_deltas = []
# for i, rand_delta in enumerate(rand_deltas):
#     # print(rand_delta)
#     for j, rand_vel in enumerate(rand_vels): 
#         solved_rand_delta = system.interact_ode(t_span, rand_delta, rand_vel, v_eqm=True, dv=True)
#         solved_rand_deltas.append(solved_rand_delta)
#         np.savetxt(f'data/rand_delta_{i}_{j}.txt', solved_rand_delta)

for i, dy in enumerate(dys):
    # print(delta)
    single_delta = system.interact_ode(t_span, dy, [0,0], v_eqm=True, dv=False)
    np.savetxt(f'data/single_delta{i}.txt', single_delta)



# measured_period = measure_period(ODE(unperturbed))
# print(f"calculated time period is {T}")
# system.planar_plot(-7, 7, -5, 9)