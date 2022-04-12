import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import math
import random

G = 4*math.pi**2
N = 10000 # number of evaluated points on the orbit
t_span = [0,50]
t_eval = np.linspace(t_span[0], t_span[1], N)
asteroids = 100

# formulaic deltas starting at (-0.001,0.001) and ending at (-0.1, 0.1)
delta_xs = np.linspace(-0.001, -0.1, asteroids, endpoint = True)
delta_ys = np.linspace(0.001, 0.1, asteroids, endpoint = True)
deltas = np.array(list(zip(delta_xs, delta_ys)))
deltas = np.array(deltas)

# randomised deltas in the square with corners (-0.1, -0.1) and (0.1, 0.1)
rand_deltas = []
for i in range(asteroids):
    rand_delta_x = (1/10000) * random.randint(-1000, 1000)
    rand_delta_y = (1/10000) * random.randint(-1000, 1000)
    rand_deltas.append([rand_delta_x, rand_delta_y])
rand_deltas = np.array(rand_deltas)





class Point:
    def __init__(self, r):
        self.x = r[0]
        self.y = r[1]


class ODE:
    def __init__(self, ode):
        self.ode = ode

        self.r_j = Point([ode[:, 0], ode[:, 2]])
        self.v_j = Point([ode[:, 1], ode[:, 3]])

        self.r_s = Point([ode[:, 4], ode[:, 6]])
        self.v_s = Point([ode[:, 5], ode[:, 7]])

        self.r_a = Point([ode[:, 8], ode[:, 10]])
        self.v_a = Point([ode[:, 9], ode[:, 11]])

class Conditions:
    def __init__(self, conditions):
        self.conditions = conditions

        self.r_j = Point([conditions[0], conditions[2]])
        self.v_j = Point([conditions[1], conditions[3]])

        self.r_s = Point([conditions[4], conditions[6]])
        self.v_s = Point([conditions[5], conditions[7]])

        self.r_a = Point([conditions[8], conditions[10]])
        self.v_a = Point([conditions[9], conditions[11]])


def format_conditions(r_j, v_j, r_s, v_s, r_a, v_a):
    '''
    puts readable data into necessary form for ODE solver
    '''
    return np.array([

        r_j[0], v_j[0], r_j[1], v_j[1],

        r_s[0], v_s[0], r_s[1], v_s[1],

        r_a[0], v_a[0], r_a[1], v_a[1],

    ])


def time_period(R, m_j, m_s):
    '''
    returns time period of an asteroid in equilibrium position 
    (same T as Jupiter around Sun)
    '''
    T = ( (4 * math.pi**2 * R**3 )/ (G * (m_s + m_j)) ) ** (1/2)
    # print(T)
    return T

def circ_orbit_conditions(y_j, m_j, m_s, delta):

    '''
    fix velocity of J and S such that they make circular orbits
    y_j fixes the initial positions of the Sun and Jupiter, as well as their initial velocities
    origin is calculated to be the CoM of the system
    delta are vectors which separate the asteroids from their equilibrium position (Trojan/Greek)

    WE ARE WORKING IN THE COM FRAME IN XY PLANE ONLY
    
    returns in the form required to enter into ode solver below
    '''

    y_s = - y_j * (m_j/m_s) # fixed such that CoM lies at origin

    r_j = [0, y_j] # initially, x is zero for both Sun and Jupiter
    r_s = [0, y_s]

    x_a = (3**(1/2) / 2) * (abs(y_s) + abs(y_j))
    y_a =  (1/2) * (abs(y_j) - abs(y_s))
    
    r_eqm = np.array([x_a, y_a])


    r_a = r_eqm + delta


    # get initial velocities of Jupiter, Sun, Asteroid at Lagrange Point

    R = abs(r_j[1]) + abs(r_s[1]) # Distance between Sun and Jupiter

    vx_s = - (G*m_s/R)**(1/2) * m_j/(m_j + m_s)
    vy_s = 0

    vx_j = (G*m_s/R)**(1/2) * m_s/(m_j + m_s) 
    vy_j = 0

    vx_a = (m_s - m_j )/ (2*(m_s + m_j)) * (G*m_s/R)**(1/2)
    vy_a = -(3*G*m_s / (4*R))**(1/2)

    v_s = [vx_s, vy_s]
    v_j = [vx_j, vy_j]
    v_a = [vx_a, vy_a]


    conditions = format_conditions(r_j, v_j, r_s, v_s, r_a, v_a)

    return conditions

def circ_orbit_conditions_new(R, m_j, m_s, delta):

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
    v_a = np.array([v * np.cos(th) , - v * np.sin(th)])


    conditions = format_conditions(r_j, v_j, r_s, v_s, r_a, v_a)

    return conditions

def circ_orbit_conditions_new_new(R, m_j, m_s, delta):

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
    v_a = np.array([v * np.cos(th) , - v * np.sin(th)])


    conditions = format_conditions(r_j, v_j, r_s, v_s, r_a, v_a)

    return conditions


def measure_period(ode):
    '''
    T/2 is when the x coordinate of Jupiter becomes negative
    '''
    index = np.where(ode.r_j.x < 0)[0][0] # first index that contains negative x value
    
    half_period = t_eval[index]
    period = 2 * half_period
    # print(f"measured period is {period}")
    return period


class Two_Body_System:

    def __init__(self, m_j, m_s):
        self.m_j = m_j 
        self.m_s = m_s 


    def general_orbit_conditions(self, r_j, v_j, r_s, v_s, r_a, v_a):
        '''
        test solver for any general orbit, not necessarily circular
        '''

        conditions = format_conditions(r_j, v_j, r_s, v_s, r_a, v_a)

        return conditions
    
        
    def interact_ivp(self, t_span, delta):
        
        '''
        m, r are the mass (scalar) and position (array) of the object respectively
        Make the usual equation of motion into coupled first order ODEs
        
        '''

        y_0 = circ_orbit_conditions_new(1.01, self.m_j, self.m_s, delta)
        # y_0 = self.general_orbit_conditions([0,1], [0,0], [0,-0.1], [0,0], [0,0], [0,0])        
        
        
        sol = scipy.integrate.solve_ivp(
            self.ODE, 
            t_span, 
            y_0, 
            method='RK45',
            t_eval = t_eval
            )
        
        return sol

    def interact_ode(self, t_span, delta):
        '''
        Using a different solver - a better one!
        '''
        
        # y_0 = circ_orbit_conditions(1, self.m_j, self.m_s, delta)
        y_0 = circ_orbit_conditions_new(1.01, self.m_j, self.m_s, delta)
        
        sol = scipy.integrate.odeint(
            self.ODE, 
            y_0, 
            t_eval,
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






    # def cutoff(self, Us, cutoff):
    #     '''
    #     cuts of the potential below a certain value 
    #     '''

    #     for U in Us:
    #         for val in U:
    #             # print(val)
    #             if val < cutoff:
    #                 val = 0
    #     return Us


    def field(self, x, y):

        # return x+y
        R = 1.01
        conditions = circ_orbit_conditions_new(R, self.m_j, self.m_s, [0,0])
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


    def planar_plot(self, x_min, x_max, y_min, y_max):

        '''
        somehow pass parameters into this to go into the field function, such as R (defined in both functions...), and ode
        '''
        # conditions = circ_orbit_conditions(1, self.m_j, self.m_s, [0,0])
        R = 1.01
        conditions = circ_orbit_conditions_new(R, self.m_j, self.m_s, [0,0])
        ode = Conditions(conditions)

        x = np.linspace(x_min, x_max, 256)
        y = np.linspace(y_min, y_max, 256)
        X, Y = np.meshgrid(x, y)
        U_eff = self.field(X, Y)

        fig, ax = plt.subplots()

        levels = np.arange(-100,-50,1);
        plt.contour(X, Y, U_eff, levels);
        ax.set_aspect('equal', adjustable='box');
        plt.plot(ode.r_a.x, ode.r_a.y , marker="x", color='k');
        plt.plot(-ode.r_a.x, ode.r_a.y , marker="x", color='k', label = "L4/L5");
        plt.plot(ode.r_j.x, ode.r_j.y , marker="o", color='b', label = 'Jupiter');
        plt.plot(ode.r_s.x, ode.r_s.y , marker="o", color='r', label = 'Sun');
        ax.set_xlabel("x");
        ax.set_ylabel("y");
        ax.legend(loc='lower right');
        plt.title("Effective Potential for co-rotating frame");
        plt.show();





# trojan_asteroid = Asteroid([0,0])
planet = Two_Body_System(0.01, 1)
unperturbed = planet.interact_ode(t_span, [0,0])
np.savetxt('unperturbed.txt', unperturbed)

solved_deltas = []
for i, delta in enumerate(deltas):
    # print(delta)
    solved_delta = planet.interact_ode(t_span, delta)
    solved_deltas.append(solved_delta)
    np.savetxt(f'delta_{i}.txt', solved_delta)


solved_rand_deltas = []
for i, rand_delta in enumerate(rand_deltas):
    # print(rand_delta)
    solved_rand_delta = planet.interact_ode(t_span, rand_delta)
    solved_rand_deltas.append(solved_rand_delta)
    np.savetxt(f'rand_delta_{i}.txt', solved_rand_delta)


# y_j = planet.circ_orbit_conditions(1, [0,0])[2]
# T = planet.time_period(unpert_ode.r_j.y[0])


# sol = planet.interact_ivp(t_span, delta).y;

measure_period(ODE(unperturbed))






        
        

    
        
