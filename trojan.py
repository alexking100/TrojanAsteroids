import numpy as np
import scipy
from scipy.integrate import solve_ivp
import math


def format_conditions(r_j, v_j, r_s, v_s, r_T, v_T):
    # puts readable data into necessary form for solver
    return[

        r_j[0], v_j[0], r_j[1], v_j[1],

        r_s[0], v_s[0], r_s[1], v_s[1],

        r_T[0], v_T[0], r_T[1], v_T[1],

    ]


class MassiveBody:
    
    
    def initial_conditions(self, y_j, delta_T, m_j, m_s):

        '''
        fix velocity of J and S such that they make circular orbits
        y_j fixes the initial positions of the Sun and Jupiter, as well as their initial velocities
        origin is calculated to be the CoM of the system
        delta_T/G are vectors which separate the asteroids from their equilibrium position (Trojan/Greek)

        WE ARE WORKING IN THE COM FRAME IN XY PLANE ONLY

        '''

        G = 4 * math.pi**2 # in the units here, G = 4*pi^2
        y_s = - y_j * (m_j/m_s) # fixed such that CoM lies at origin

        r_j = [0, y_j] # initially, x is zero for both Sun and Jupiter
        r_s = [0, y_s]


        # eqm_trojan is the equilibrium position of a Trojan asteroid
        # eqm_greek is the equilibrium position of a Greek asteroid

        trojan_x = (3**(1/2) / 2) * (abs(r_s[1]) + abs(r_j[1]))
        trojan_y =  (1/2) * (abs(r_j[1]) - abs(r_s[1]))
        # greek_x = - trojan_x
        # greek_y = trojan_y
        
        eqm_trojan = [trojan_x, trojan_y]
        # eqm_greek = [greek_x, greek_y]


        r_T = eqm_trojan + delta_T
        # r_G = eqm_greek + delta_G
        # HOLD OFF ON GREEK ASTEROIDS FOR NOW...


        R = abs(r_j[1]) + abs(r_s[1]) # Distance between Sun and Jupiter

        vx_s = - (G*m_s/R)**(1/2) * m_j/(m_j + m_s)
        vy_s = 0

        vx_j = (G*m_s/R)**(1/2) * m_s/(m_j + m_s)
        vy_j = 0

        vx_T = (m_s - m_j )/ (2*(m_s + m_j)) * (G*m_s/R)**(1/2)
        vy_T = -(3*G*m_s / (4*R))**(1/2)

        v_s = [vx_s, vy_s]
        v_j = [vx_j, vy_j]
        v_T = [vx_T, vy_T]

        # v_T = [(G*m_s/R)**(1/2) * 1/2, - (3*G*m_s/R)**(1/2) * 1/2]
        # v_G = [(m_s - m_j )/ (2*(m_s + m_j)) * (G*m_s/R), (3*G*m_s / (4*R))**(1/2)]

        conditions = format_conditions(r_j, v_j, r_s, v_s, r_T, v_T)

        return conditions

        
    def interact(self, m_j, m_s, t_span):
        
        '''
        m, r are the mass (scalar) and position (array) of the object respectively
        Make the usual equation of motion into coupled first order ODEs
        
        '''

        y_0 = self.initial_conditions(1, [0,0], m_j, m_s)

        params = [m_j, m_s]
        t_eval = np.linspace(t_span[0], t_span[1], 10000)
        
        
        
        sol = scipy.integrate.solve_ivp(
            self.ODE, 
            t_span, 
            y_0, 
            method='RK45',
            t_eval = t_eval,
            args=params
            )
        
        return sol
    
    
    def ODE(self, t, y, m_j, m_s):
        
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
            G * m_s * dx_js / (dx_js**2 + dy_js**2)**(3/2),
            G * m_s * dy_js / (dx_js**2 + dy_js**2)**(3/2)
        ]

        accel_s = [
            - G * m_j * dx_js / (dx_js**2 + dy_js**2)**(3/2),
            - G * m_j * dy_js / (dx_js**2 + dy_js**2)**(3/2)
        ]

        accel_T_x = G * m_j * dx_Tj / (dx_Tj**2 + dy_Tj**2)**(3/2) + G * m_s * dx_Ts / (dx_Ts**2 + dy_Ts**2)**(3/2)
        accel_T_y = G * m_j * dy_Tj / (dx_Tj**2 + dy_Tj**2)**(3/2) + G * m_s * dy_Ts / (dx_Ts**2 + dy_Ts**2)**(3/2)

        accel_T = [accel_T_x, accel_T_y]
        
        conditions = format_conditions(v_j, accel_j, v_s, accel_s, v_T, accel_T)

        return conditions
    
# class Asteroid():

#     def __init__(self, delta, trojan = True):
#         '''
#         delta_T/G are initial positions of asteroids with respect to Trojan / Greek equilibrium positions
#         '''
#         self.delta = delta
#         self.trojan = trojan

#     def positions(r_s, r_j):

#         # eqm_trojan is the equilibrium position of a Trojan asteroid
#         # eqm_greek is the equilibrium position of a Greek asteroid

#         if trojan:
#             eqm_trojan = [(3**(1/2) / 2) * (abs(r_s[1]) + abs(r_j[1])), -(1/2) * (abs(r_s[1]) - abs(r_j[1]))]
#             r = eqm_trojan + self.delta
#         else:
#             eqm_greek = [- (3**(1/2) / 2) * (abs(r_s[1]) + abs(r_j[1])), -(1/2) * (abs(r_s[1]) - abs(r_j[1]))]
#             r = eqm_greek + self.delta

#         return r







# trojan_asteroid = Asteroid([0,0])
planet = MassiveBody()

# set parameters of system
m_j = 0.001
m_s = 1
t_span = [0,10]


sol = planet.interact(m_j, m_s, t_span);

np.savetxt('orbitdata.txt', sol.y)


        
        

    
        
