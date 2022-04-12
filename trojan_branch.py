import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import math

G = 4*math.pi**2

def format_conditions(r_j, v_j, r_s, v_s, r_T, v_T):
    # puts readable data into necessary form for solver
    return[

        r_j[0], v_j[0], r_j[1], v_j[1],

        r_s[0], v_s[0], r_s[1], v_s[1],

        r_T[0], v_T[0], r_T[1], v_T[1],

    ]


def r_asteroid(r_j, r_s, delta = np.zeros(2), trojan = True):

    # eqm_trojan is the equilibrium position of a Trojan asteroid
    # eqm_greek is the equilibrium position of a Greek asteroid

    eqm_trojan_x = (3**(1/2) / 2) * (abs(r_s[1]) + abs(r_j[1]))
    eqm_trojan_y = -(1/2) * (abs(r_s[1]) - abs(r_j[1]))

    eqm_greek_x = - (3**(1/2) / 2) * (abs(r_s[1]) + abs(r_j[1]))
    eqm_greek_y = - (1/2) * (abs(r_s[1]) - abs(r_j[1]))

    if trojan:
        eqm_trojan = np.array([eqm_trojan_x, eqm_trojan_y])
        r = eqm_trojan + self.delta
    else:
        eqm_greek = np.array([eqm_greek_x, eqm_greek_y])
        r = eqm_greek + self.delta

    return r

def r_accel()




def v_asteroid(r_j, r_s, m_j, m_s, delta, trojan = True):

    R = abs(r_j[1]) + abs(r_s[1]) # Distance between Sun and Jupiter

    if trojan:
        vx = (m_s - m_j )/ (2*(m_s + m_j)) * (G*m_s/R)**(1/2)
        vy = -(3*G*m_s / (4*R))**(1/2)

        v = np.array([vx, vy])

    return v


class Two_Body_System:

    def __init__(self, m_j, m_s):
        self.m_j = m_j 
        self.m_s = m_s

    def general_orbit_conditions(self, r_j, v_j, r_s, v_s, r_T, v_T):
        '''
        test solver for any general orbit, not necessarily circular
        '''

        conditions = format_conditions(r_j, v_j, r_s, v_s, r_T, v_T)

        return conditions
    
    
    def circ_orbit_conditions(self, y_j, delta = np.zeros(2)):

        '''
        fix velocity of J and S such that they make circular orbits
        y_j fixes the initial positions of the Sun and Jupiter, as well as their initial velocities
        origin is calculated to be the CoM of the system
        Asteroid.delta are vectors which separate the asteroids from their equilibrium position (Trojan/Greek)

        WE ARE WORKING IN THE COM FRAME IN XY PLANE ONLY

        '''

        y_s = - y_j * (self.m_j/self.m_s) # fixed such that CoM lies at origin

        r_j = [0, y_j] # initially, x is zero for both Sun and Jupiter
        r_s = [0, y_s]

        r_asteroid = r_asteroid(r_j, r_s, delta)

        print(delta)
        print(r_asteroid)



        R = abs(r_j[1]) + abs(r_s[1]) # Distance between Sun and Jupiter

        vx_s = - (G*self.m_s/R)**(1/2) * self.m_j/(self.m_j + self.m_s)
        vy_s = 0

        vx_j = (G*self.m_s/R)**(1/2) * self.m_s/(self.m_j + self.m_s)
        vy_j = 0

        v_s = [vx_s, vy_s]
        v_j = [vx_j, vy_j]

        v_asteroid = v_asteroid(r_j, r_s, self.m_j, self.m_s, delta)

        conditions = format_conditions(r_j, v_j, r_s, v_s, r_asteroid, v_asteroid)

        return conditions

        
    def interact_ivp(self, t_span, delta):
        
        '''
        m, r are the mass (scalar) and position (array) of the object respectively
        Make the usual equation of motion into coupled first order ODEs
        
        '''

        y_0 = self.circ_orbit_conditions(1, delta)
        # y_0 = self.general_orbit_conditions([0,1], [0,0], [0,-0.1], [0,0], [0,0], [0,0])

        t_eval = np.linspace(t_span[0], t_span[1], 10000)
        
        
        
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
        Using a different solver 
        '''

        y_0 = self.circ_orbit_conditions(1, delta)
        # y_0 = self.general_orbit_conditions([0,1], [0,0], [0,-0.1], [0,0], [0,0], [0,0])

        t_eval = np.linspace(t_span[0], t_span[1], 10000)
        
        
        sol = scipy.integrate.odeint(
            self.ODE, 
            y_0, 
            t_eval,
            tfirst = True
        )
        

        return sol
    
    def ODE(self, t, y):
        
#       y = [x_j, vx_j, y_j, vy_j, x_s, vx_s, y_s, vy_s, x_a, vx_a, y_a, vy_a]

        G = 4 * math.pi**2 # in the units here, G = 4*pi^2

        # assign meaning to y-values
        r_j = [y[0], y[2]]
        v_j = [y[1], y[3]]

        r_s = [y[4], y[6]]
        v_s = [y[5], y[7]]

        r_T = [y[8], y[10]]
        v_T = [y[9], y[11]]
        
        # convention is dx_js is the x vector from Jupiter to the Sun, etc.
        dx_js = r_s[0] - r_j[0]
        dx_Ts = r_s[0] - r_T[0]
        dx_Tj = r_j[0] - r_T[0]

        dy_js = r_s[1] - r_j[1]
        dy_Ts = r_s[1] - r_T[1]
        dy_Tj = r_j[1] - r_T[1]


        accel_j = [
            G * self.m_s * dx_js / (dx_js**2 + dy_js**2)**(3/2),
            G * self.m_s * dy_js / (dx_js**2 + dy_js**2)**(3/2)
        ]

        accel_s = [
            - G * self.m_j * dx_js / (dx_js**2 + dy_js**2)**(3/2),
            - G * self.m_j * dy_js / (dx_js**2 + dy_js**2)**(3/2)
        ]

        accel_T_x = G * self.m_j * dx_Tj / (dx_Tj**2 + dy_Tj**2)**(3/2) + G * self.m_s * dx_Ts / (dx_Ts**2 + dy_Ts**2)**(3/2)
        accel_T_y = G * self.m_j * dy_Tj / (dx_Tj**2 + dy_Tj**2)**(3/2) + G * self.m_s * dy_Ts / (dx_Ts**2 + dy_Ts**2)**(3/2)

        accel_T = [accel_T_x, accel_T_y]
        
        conditions = format_conditions(v_j, accel_j, v_s, accel_s, v_T, accel_T)

        return conditions

    def deviation():
        '''
        difference between deviated asteroid and asteroid at equilibrium
        '''





    # def energy():

    #     '''
    #     calculate the total energy of a given object
    #     '''

    #     potential = -G*self.m_s











# trojan_asteroid = Asteroid([0,0])
planet = Two_Body_System(0.01, 1)

t_span = [0,10]

delta_T = [0,0]

sol = planet.interact_ode(t_span, delta_T);
# sol = planet.interact_ivp(t_span, delta_T).y;

np.savetxt('orbitdata.txt', sol)


        
        

    
        
