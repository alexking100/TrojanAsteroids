import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import math
import random

'''
Setting up parameters
'''

G = 4*math.pi**2 # value of Gravitational Constant in our units
N = 1000 # number of evaluated points on the orbit
t_span = [0,50] # time span of orbit
t = np.linspace(t_span[0], t_span[1], N) 
asteroids = 1 # number of asteroids in the orbit
m_j = 0.01 # in units of solar mass
m_s = 1 # in units of solar mass
R = 5.2 # distance from Jupiter to Sun, in AU 
max_delta = 1
min_delta = 0.1
max_delta_v = 0.1
min_delta_v = 0

'''
Define spatial deviation from lagrange point for each asteroid (formulaic and random)
'''

# formulaic deltas starting at (-0.001,0.001) and ending at (-0.1, 0.1)
# direction changes in the y direction
# velocity changes (for each y) in the x direction

def get_ds(asteroids):
    return np.linspace(-max_delta, max_delta, asteroids)


def get_deltas(asteroids):
    delta_xs = np.zeros(asteroids)
    delta_ys = np.linspace(min_delta, max_delta, asteroids, endpoint = True)
    deltas = np.array(list(zip(delta_xs, delta_ys)))
    return deltas

def get_delta_vs(asteroids):
    delta_vxs = np.linspace(min_delta, max_delta_v, asteroids, endpoint = True)
    delta_vys = np.zeros(asteroids)
    delta_vs = np.array(list(zip(delta_vxs, delta_vys)))
    return delta_vs

# Randomised deltas in the square with corners (-0.1, -0.1) and (0.1, 0.1)
# added to eqm positions to test stability of orbit
def get_rand_deltas(asteroids):
    rand_deltas = []
    for i in range(asteroids):
        rand_delta_x = max_delta/1000 * random.randint(-1000, 1000)
        rand_delta_y = min_delta/1000 * random.randint(-1000, 1000)
        rand_deltas.append([rand_delta_x, rand_delta_y])
    rand_deltas = np.array(rand_deltas)
    return rand_deltas

# Random (small) velocities, added to the eqm velocity to test stability of orbit
def get_rand_vels(asteroids):
    rand_vels = []
    for i in range(asteroids):
        rand_v_x = (1/10000) * random.randint(-1000,1000)
        rand_v_y = (1/10000) * random.randint(-1000,1000)
        rand_vels.append([rand_v_x, rand_v_y])
    rand_vels = np.array(rand_vels)
    return rand_vels


# def circ_orbit_conditions(y_j, m_j, m_s, delta):

#     '''
#     fix velocity of J and S such that they make circular orbits
#     y_j fixes the initial positions of the Sun and Jupiter, as well as their initial velocities
#     origin is calculated to be the CoM of the system
#     delta are vectors which separate the asteroids from their equilibrium position (Trojan/Greek)

#     WE ARE WORKING IN THE COM FRAME IN XY PLANE ONLY
    
#     returns in the form required to enter into ode solver below
#     '''

#     y_s = - y_j * (m_j/m_s) # fixed such that CoM lies at origin

#     r_j = [0, y_j] # initially, x is zero for both Sun and Jupiter
#     r_s = [0, y_s]

#     x_a = (3**(1/2) / 2) * (abs(y_s) + abs(y_j))
#     y_a =  (1/2) * (abs(y_j) - abs(y_s))
    
#     r_eqm = np.array([x_a, y_a])


#     r_a = r_eqm + delta


#     # get initial velocities of Jupiter, Sun, Asteroid at Lagrange Point

#     R = abs(r_j[1]) + abs(r_s[1]) # Distance between Sun and Jupiter

#     vx_s = - (G*m_s/R)**(1/2) * m_j/(m_j + m_s)
#     vy_s = 0

#     vx_j = (G*m_s/R)**(1/2) * m_s/(m_j + m_s) 
#     vy_j = 0

#     vx_a = (m_s - m_j )/ (2*(m_s + m_j)) * (G*m_s/R)**(1/2)
#     vy_a = -(3*G*m_s / (4*R))**(1/2)

#     v_s = [vx_s, vy_s]
#     v_j = [vx_j, vy_j]
#     v_a = [vx_a, vy_a]


#     conditions = format_conditions(r_j, v_j, r_s, v_s, r_a, v_a)

#     return conditions


def circ_orbit_conditions_new(R, m_j, m_s, delta, rand_vel, eqm = False):
    '''
    returns initial conditions for all objects (J, S, Asteroid)
    to be in circular orbit about the CoM (which lies at the origin)
    if eqm is set to True, the function calculates the correct asteroid velocity to maintain circ orbit
    if eqm is set to False, the function takes the input rand_vel as the asteroid's initial velocity
    '''

    M = m_j + m_s # total mass
    gamma = m_j / M

    r_j = np.array([0, (1 - gamma) * R])
    r_s = np.array([0, - gamma * R])

    mu = (G * M / R) **(1/2)

    v_j = np.array([(1-gamma) * mu, 0])
    v_s = np.array([- gamma * mu, 0])


    X = R * (1 - gamma + gamma**2)**(1/2) # distance between origin and Lagrange point
    th = (1 + gamma) * math.pi / 3 # angle between line to Lagrange point and y axis

    T = time_period(R, m_j, m_s)

    # print(f"calculated time period is {T}")
    omega = 2*math.pi / T
    v = omega * X

    # asteroid conditions
    r_a_eqm = np.array([X * np.sin(th), X * np.cos(th)])
    r_a = r_a_eqm + delta

    if eqm:
        v_a = np.array([v * np.cos(th) , - v * np.sin(th)])
    else:
        v_a = np.array([v * np.cos(th) , - v * np.sin(th)]) + np.array(rand_vel)

    conditions = format_conditions(r_j, v_j, r_s, v_s, r_a, v_a)

    return conditions


def time_period(R, m_j, m_s):
    '''
    returns time period of an asteroid in equilibrium position 
    (same T as Jupiter around Sun)
    '''
    T = ( (4 * math.pi**2 * R**3 )/ (G * (m_s + m_j)) ) ** (1/2)
    # print(T)
    return T


def format_conditions(r_j, v_j, r_s, v_s, r_a, v_a):
    '''
    puts readable data into necessary form for ODE solver
    '''
    return np.array([

        r_j[0], v_j[0], r_j[1], v_j[1],

        r_s[0], v_s[0], r_s[1], v_s[1],

        r_a[0], v_a[0], r_a[1], v_a[1],

    ])


class Point:
    '''
    Abstract Class
    '''
    def __init__(self, r):
        self.x = r[0]
        self.y = r[1]

class ODE:
    '''
    Abstract Class
    input is solution to ode solver
    '''
    def __init__(self, ode):
        self.ode = ode

        self.r_j = Point([ode[:, 0], ode[:, 2]])
        self.v_j = Point([ode[:, 1], ode[:, 3]])

        self.r_s = Point([ode[:, 4], ode[:, 6]])
        self.v_s = Point([ode[:, 5], ode[:, 7]])

        self.r_a = Point([ode[:, 8], ode[:, 10]])
        self.v_a = Point([ode[:, 9], ode[:, 11]])

class Conditions:
    '''
    Abstract Class
    input is array in form returned by format_conditions()
    '''
    def __init__(self, conditions):
        self.conditions = conditions

        self.r_j = Point([conditions[0], conditions[2]])
        self.v_j = Point([conditions[1], conditions[3]])

        self.r_s = Point([conditions[4], conditions[6]])
        self.v_s = Point([conditions[5], conditions[7]])

        self.r_a = Point([conditions[8], conditions[10]])
        self.v_a = Point([conditions[9], conditions[11]])


class Two_Body_System:
    '''
    creates an orbit of Jupiter/Sun/Asteroid around the CoM of the system
    can solve the differential equations
    can find and plot the effective potential about the initial setup
    '''

    def __init__(self, m_j, m_s):
        self.m_j = m_j 
        self.m_s = m_s 


    def general_orbit_conditions(self, r_j, v_j, r_s, v_s, r_a, v_a):
        '''
        test solver for any general orbit, not necessarily circular
        gives complete control for me to vary system as needed for testing
        '''

        conditions = format_conditions(r_j, v_j, r_s, v_s, r_a, v_a)

        return conditions
    
        
    def interact_ivp(self, t_span, t, method, delta, rand_vel, eqm = False):
        '''
        m, r are the mass (scalar) and position (array) of the object respectively
        Make the usual equation of motion into coupled first order ODEs
        '''

        y_0 = circ_orbit_conditions_new(R, self.m_j, self.m_s, delta, rand_vel, eqm = eqm)
        # y_0 = self.general_orbit_conditions([0,1], [0,0], [0,-0.1], [0,0], [0,0], [0,0])        
        
        
        sol = scipy.integrate.solve_ivp(
            self.ODE, 
            t_span, 
            y_0, 
            method=method,
            t_eval = t,
            )
        
        return sol


    def interact_ode(self, t_span, delta, rand_vel, eqm = False):
        '''
        Using a different solver
        '''
        
        # y_0 = circ_orbit_conditions(1, self.m_j, self.m_s, delta)
        y_0 = circ_orbit_conditions_new(R, self.m_j, self.m_s, delta, rand_vel, eqm = eqm)
        
        sol = scipy.integrate.odeint(
            self.ODE, 
            y_0, 
            t,
            tfirst = True
        )

        return sol

    
    def ODE(self, t, y):
        
#       y = [x_j, vx_j, y_j, vy_j, x_s, vx_s, y_s, vy_s, x_a, vx_a, y_a, vy_a]

        # assign meaning to y-values
        r_j = np.array([y[0], y[2]])
        v_j = np.array([y[1], y[3]])

        r_s = np.array([y[4], y[6]])
        v_s = np.array([y[5], y[7]])

        r_a = np.array([y[8], y[10]])
        v_a = np.array([y[9], y[11]])
        
        # convention is dx_js is the x vector from Jupiter to the Sun, etc.
        dx_js = r_s[0] - r_j[0]
        dx_as = r_s[0] - r_a[0]
        dx_aj = r_j[0] - r_a[0]

        dy_js = r_s[1] - r_j[1]
        dy_as = r_s[1] - r_a[1]
        dy_aj = r_j[1] - r_a[1]


        accel_j = np.array([
            G * self.m_s * dx_js / (dx_js**2 + dy_js**2)**(3/2),
            G * self.m_s * dy_js / (dx_js**2 + dy_js**2)**(3/2)
        ])

        accel_s = np.array([
            - G * self.m_j * dx_js / (dx_js**2 + dy_js**2)**(3/2),
            - G * self.m_j * dy_js / (dx_js**2 + dy_js**2)**(3/2)
        ])

        accel_a_x = G * self.m_j * dx_aj / (dx_aj**2 + dy_aj**2)**(3/2) + G * self.m_s * dx_as / (dx_as**2 + dy_as**2)**(3/2)
        accel_a_y = G * self.m_j * dy_aj / (dx_aj**2 + dy_aj**2)**(3/2) + G * self.m_s * dy_as / (dx_as**2 + dy_as**2)**(3/2)

        accel_a = np.array([accel_a_x, accel_a_y])
        
        conditions = format_conditions(v_j, accel_j, v_s, accel_s, v_a, accel_a)

        return conditions


    def field(self, x, y):

        # return x+y
        conditions = circ_orbit_conditions_new(R, self.m_j, self.m_s, [0,0], [0,0], eqm = True)
        T = time_period(R, m_j, m_s)
        omega = 2 * math.pi / T
        ode = Conditions(conditions)

        r_j_array = np.array([ode.r_j.x, ode.r_j.y])
        r_s_array = np.array([ode.r_s.x, ode.r_s.y])

        r = (x**2 + y**2)**(1/2)
        delta_x_j = (x - ode.r_j.x)
        delta_y_j = (y - ode.r_j.y)
        delta_r_j = (delta_x_j**2 + delta_y_j**2)**(1/2)

        delta_x_s = (x - ode.r_s.x)
        delta_y_s = (y - ode.r_s.y)
        delta_r_s = (delta_x_s**2 + delta_y_s**2)**(1/2)

        U_j = - G * self.m_j / delta_r_j
        U_s = - G * self.m_s / delta_r_s
        U_rot = - 1/2 * r**2 * omega**2

        U_eff = np.array(U_rot) + np.array(U_j) + np.array(U_s)

        # U_eff = self.cutoff(U_eff, -50)

        return U_eff


    def planar_plot(self, x_min, x_max, y_min, y_max, solved_deltas, unpert_ode):

        '''
        somehow pass parameters into this to go into the field function, such as R (defined in both functions...), and ode
        '''
        # conditions = circ_orbit_conditions(1, self.m_j, self.m_s, [0,0])
        conditions = circ_orbit_conditions_new(R, self.m_j, self.m_s, [0,0], [0,0], eqm = True)
        ode = Conditions(conditions)

        x = np.linspace(x_min, x_max, 256)
        y = np.linspace(y_min, y_max, 256)
        X, Y = np.meshgrid(x, y)
        U_eff = self.field(X, Y)

        fig, ax = plt.subplots()

        levels = np.arange(-30, 0, 0.5); # plot min U_eff, max U_eff, step in U_eff
        ax.contour(X, Y, U_eff, levels); #plt? 
        ax.set_aspect('equal', adjustable='box');
        
        # for i, solved_delta in enumerate(solved_deltas):
        #     r_primed_delta = perform_rotation_ode(unpert_ode, solved_delta)
        #     initial_x = r_primed_delta[0][0]
        #     initial_y = r_primed_delta[1][0]
        #     if i == 0:
        #         ax.plot(initial_x, initial_y, marker = 'o', color = 'r', label = "Displaced Asteroid")
        #     else:
        #         ax.plot(initial_x, initial_y, marker = 'o', color = 'r')

        
        ax.plot(ode.r_a.x, ode.r_a.y , marker="x", color='k'); 
        ax.plot(-ode.r_a.x, ode.r_a.y , marker="x", color='k', label = "L4/L5"); 
        ax.plot(ode.r_j.x, ode.r_j.y , marker="o", markersize=10, color='b', label = 'Jupiter'); 
        ax.plot(ode.r_s.x, ode.r_s.y , marker="o", markersize=15, color='y', label = 'Sun'); 
        ax.set_xlabel("x /AU");
        ax.set_ylabel("y /AU");
        ax.legend(loc='lower right');
        ax.set_title("Effective Potential for Co-Rotating Frame"); 



def measure_period(ode):
    '''
    T/2 is when the x coordinate of Jupiter becomes negative
    To check reliability of calculated time_period() function
    '''
    index = np.where(ode.r_j.x < 0)[0][0] # first index that contains negative x value
    
    half_period = t[index]
    period = 2 * half_period
    # print(f"measured period is {period}")
    return period


def make_plots_ode(unpert_ode, solved_deltas):

    fig, ax = plt.subplots()
    ax.plot(unpert_ode.r_j.x, unpert_ode.r_j.y, label = "Jupiter")
    ax.plot(unpert_ode.r_s.x, unpert_ode.r_s.y, label = "Sun")
    ax.plot(unpert_ode.r_a.x, unpert_ode.r_a.y, label = "Unperturbed Asteroid")
    ax.set_xlim(-8,8)
    ax.set_ylim(-8,8)
    ax.plot(solved_deltas[0].r_a.x, solved_deltas[0].r_a.y, label = "Perturbed Asteroid")
    ax.legend()
    ax.set_title("Motion of Sun, Jupiter and Asteroids in non-rotating frame")


def stationary_initial_position(unpert_ode, solved_deltas):

    r_primed_unpert = perform_rotation_ode(unpert_ode, unpert_ode)
    unpert_x = r_primed_unpert[0][0]
    unpert_y = r_primed_unpert[1][0]

    fig, ax = plt.subplots()
    for i, solved_delta in enumerate(solved_deltas):
        r_primed_delta = perform_rotation_ode(unpert_ode, solved_delta)
        initial_x = r_primed_delta[0][0]
        initial_y = r_primed_delta[1][0]
        ax.plot(initial_x, initial_y, marker = 'o', color = 'r', label = "Displaced Asteroid")

    ax.plot(unpert_x, unpert_y, marker = 'o', color = 'b', label = "Asteroid at L4")
    ax.plot(unpert_ode.r_j.x[0], unpert_ode.r_j.y[0], marker = 'x', label = "Jupiter")
    ax.plot(unpert_ode.r_s.x[0], unpert_ode.r_s.y[0], marker = 'x', label = "Sun")
    ax.set_xlim(-4, 4)
    ax.set_ylim(-2, 6)
    ax.legend()
    ax.set_title("Initial Positions of Bodies")


def angle(unpert_ode, t):

    T = measure_period(unpert_ode)

    omega = 2*math.pi / T
    thetas = np.zeros(len(t))

    for i, t in enumerate(t):
        thetas[i] = omega*t

    return thetas



def perform_rotation_ode(unpert_ode, ode):
    '''
    rotate the x-y frame at the same rate as the eqm asteroid
    returns two vectors in a tuple, each length N and containing x/y coordinates in the co-rotated frame
    '''
    t = np.linspace(t_span[0], t_span[1], N)
    thetas = angle(unpert_ode, t)

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


def corotated_deviation(unpert_ode, solved_deltas, deltas, rand_vels, title):

    print(len(solved_deltas))

    r_primed_unpert = perform_rotation_ode(unpert_ode, unpert_ode)

    fig, ax = plt.subplots()
    # ax.plot(r_primed_unpert[0], r_primed_unpert[1], label = "Undisplaced Asteroid")
    ax.plot(unpert_ode.r_j.x[0], unpert_ode.r_j.y[0], marker = 'o', ms = 10)
    ax.plot(unpert_ode.r_s.x[0], unpert_ode.r_s.y[0], marker = 'o', ms = 15, color = 'y')
    ax.plot(unpert_ode.r_a.x[0], unpert_ode.r_a.y[0], marker = 'x', ms = 15, color = 'k')
    # asteroid_orbit = plt.Circle((0, 0), (unpert_ode.r_a.x**2 + unpert_ode.r_a.y**2)**0.5, color='k', fill=False)
    # ax.add_patch(asteroid_orbit) this doesn't work...

    for i, solved_delta in enumerate(solved_deltas):
        r_primed_delta = perform_rotation_ode(unpert_ode, solved_delta)
        # diff_xs[i] = r_primed_delta[0] - r_primed_unpert[0]
        # diff_ys[i] = r_primed_delta[1] - r_primed_unpert[1]
        ax.plot(r_primed_delta[0], r_primed_delta[1], label = "Displaced Asteroid") #, label = f"delta_r = {float(round(np.linalg.norm(deltas[i]), 3))}" , delta_v = {float(round(np.linalg.norm(rand_vels[i]), 3))}")
    ax.set_xlim(-6,6)
    ax.set_ylim(-3,8)
    ax.legend()
    ax.set_title(title)


    # diff_ast_x = r_primed_unpert[0] - unpert_ode.r_a.x[0]
    # diff_ast_y = r_primed_unpert[1] - unpert_ode.r_a.y[0]
    # print(unpert_ode.r_a.x[0], unpert_ode.r_a.y[0])

    # plt.plot(unpert_ode.r_a.x, unpert_ode.r_a.y, label = "Asteroid in rest frame")
    # plt.plot(r_primed_unpert[0], r_primed_unpert[1], label = "Asteroid in rotated frame")
    # plt.legend()

    # fig, ax = plt.subplots()
    # plt.plot(diff_x, diff_y, label = 'Deviation of asteroid') # difference between the asteroid's position and the initial lagrange point

def energy(ode, m_j, m_s):

    v_j = (ode.v_j.x**2 + ode.v_j.y**2)**0.5
    v_s = (ode.v_s.x**2 + ode.v_s.y**2)**0.5

    r_js = ((ode.r_s.x - ode.r_j.x)**2 + (ode.r_s.y - ode.r_j.y)**2)**0.5

    print(f"vel jupiter = {v_j}")
    print(f"vel sun = {v_s}")
    print(f"sun to jupiter = {r_js}")

    T_j = (1/2) * m_j * v_j**2
    U_j = - (G * m_s * m_j) / r_js

    T_s = (1/2) * m_s * v_s**2
    U_s = - (G * m_s * m_j) / r_js

    E = T_j + U_j + T_s + U_s
    print(E)

    return E


def maximum_deviation(solved_rand_deltas, rand_deltas, unpert_ode):

    max_mod_delta_rs = []

    for r in solved_rand_deltas:
        r_primed = perform_rotation_ode(unpert_ode, r)
        r_primed_unpert = perform_rotation_ode(unpert_ode, unpert_ode)
        delta_r = r_primed - r_primed_unpert
        # print(delta_r)
        mod_delta_r = np.linalg.norm(delta_r, axis = 0)
        # print(mod_delta_r)
        # print(f"shape is {np.shape(mod_delta_r)}")
        max_mod_delta_r = np.amax(mod_delta_r)
        max_mod_delta_rs.append(max_mod_delta_r)
    mod_deltas = np.linalg.norm(rand_deltas, axis = 1)
    mod_deltas = np.repeat(mod_deltas, asteroids) 
    '''
    np.repeat function turns the array from (0,1,2,3) into (0,0,0,1,1,1,2,2,2,3,3,3) where each element
    in the array is repeated (asteroids) times.  This is necessary because there are 10 velocities for each mod_delta value.
    '''

    # print(max_mod_delta_rs)
    # print(mod_deltas)

    fig, ax = plt.subplots()
    plt.plot(mod_deltas, max_mod_delta_rs, marker = 'o')
    plt.title("Max deviation vs initial displacement from L4")
    # ax.set_ylim(0, 10)
    # return (mod_deltas, max_mod_delta_rs)


def heatmap(sols, unpert_ode, xlabel, ylabel):

    max_mod_delta_rs = []
    arr = np.zeros(asteroids**2)
    i = 0
    for ode in sols:
        r_primed = perform_rotation_ode(unpert_ode, ode)
        r_primed_unpert = perform_rotation_ode(unpert_ode, unpert_ode)
        delta_r = r_primed - r_primed_unpert
        mod_delta_r = np.linalg.norm(delta_r, axis = 0)
        max_mod_delta_r = np.amax(mod_delta_r)
        max_mod_delta_rs.append(max_mod_delta_r)

        arr[i] = max_mod_delta_r
        i += 1

    X = np.reshape(arr, (asteroids, asteroids))

    fig, ax = plt.subplots()
    ax.imshow(X, vmin = 0, vmax = 5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Maximum Deviation from L4")




        
        

    
        
