from trojan import *

deltas = get_deltas(asteroids)
delta_vs = get_delta_vs(asteroids)
rand_deltas = get_rand_deltas(asteroids)
rand_vels = get_rand_vels(asteroids)
dxs = get_ds(asteroids)
dys = get_ds(asteroids)
dv_xs = get_ds(asteroids)
dv_ys = get_ds(asteroids)

# cannot run the random functions again in another file as they will be different
# save them to .txt and call the data in another file
np.savetxt('data/rand_vels.txt', rand_vels)

np.savetxt('data/rand_deltas.txt', rand_deltas)

system = Two_Body_System(m_j, m_s)

unperturbed = system.interact_ode(t_span, [0,0], [0,0], eqm = True)
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

# get data for a grid of displaced asteroids in the x-y plane about L4
solved_grids = []
for i, dx in enumerate(dxs):
    for j, dy in enumerate(dys):
        solved_grid = system.interact_ode(t_span, [dx,dy], [0,0], eqm = False)
        solved_grids.append(solved_grid)
        np.savetxt(f'data/grid_{i}_{j}.txt', solved_grid)

solved_vs = []
for i, dv_x in enumerate(dv_xs):
    for j, dv_y in enumerate(dv_ys):
        solved_v = system.interact_ode(t_span, [0,0], [dv_x, dv_y], eqm = False)
        solved_vs.append(solved_v)
        np.savetxt(f'data/v_{i}_{j}.txt', solved_v)

# get data for formulaically displaced asteroids in y, v_x space
solved_deltas = []
for i, delta in enumerate(deltas):
    # print(delta)
    for j, delta_v in enumerate(delta_vs): 
        solved_delta = system.interact_ode(t_span, delta, delta_v, eqm = False)
        solved_deltas.append(solved_delta)
        np.savetxt(f'data/delta_{i}_{j}.txt', solved_delta)

# get data for randomly displaced asteroids (each delta value has every random velocity attached)
# solved_rand_deltas is of size (asteroids*asteroids)
solved_rand_deltas = []
for i, rand_delta in enumerate(rand_deltas):
    # print(rand_delta)
    for j, rand_vel in enumerate(rand_vels): 
        solved_rand_delta = system.interact_ode(t_span, rand_delta, rand_vel, eqm = False)
        solved_rand_deltas.append(solved_rand_delta)
        np.savetxt(f'data/rand_delta_{i}_{j}.txt', solved_rand_delta)

measured_period = measure_period(ODE(unperturbed))
calc_period = time_period(R, m_j, m_s)

# solver_methods()
# system.planar_plot(-7, 7, -5, 9)