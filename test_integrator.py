'''
script that tests the stability and accuracy of the integrators
'''

from trojan import *

unpert_ode = ODE(np.loadtxt('data/unperturbed.txt'))
unpert_ivp_RK45 = ODE(np.loadtxt('data/solver/unperturbed_ivp_RK45.txt'))
unpert_ivp_RK23 = ODE(np.loadtxt('data/solver/unperturbed_ivp_RK23.txt'))
unpert_ivp_DOP853 = ODE(np.loadtxt('data/solver/unperturbed_ivp_DOP853.txt'))
unpert_ivp_Radau = ODE(np.loadtxt('data/solver/unperturbed_ivp_Radau.txt'))
unpert_ivp_BDF = ODE(np.loadtxt('data/solver/unperturbed_ivp_BDF.txt'))
unpert_ivp_LSODA = ODE(np.loadtxt('data/solver/unperturbed_ivp_LSODA.txt'))


E_ode = energy(unpert_ode, m_j, m_s)
E_ivp_RK45 = energy(unpert_ivp_RK45, m_j, m_s)
E_ivp_RK23 = energy(unpert_ivp_RK23, m_j, m_s)
E_ivp_DOP853 = energy(unpert_ivp_DOP853, m_j, m_s)
E_ivp_Radau= energy(unpert_ivp_Radau, m_j, m_s)
E_ivp_BDF = energy(unpert_ivp_BDF, m_j, m_s)
E_ivp_LSODA = energy(unpert_ivp_LSODA, m_j, m_s)


fig, ax = plt.subplots()
ax.plot(t, E_ode, label = "odeint")
# ax.plot(t, E_ivp_RK45, label = "solve_ivp: RK45")
ax.plot(t, E_ivp_RK23, label = "solve_ivp: RK23")
ax.plot(t, E_ivp_DOP853, label = "solve_ivp: DOP853")
ax.plot(t, E_ivp_Radau, label = "solve_ivp: Radau")
ax.plot(t, E_ivp_BDF, label = "solve_ivp: BDF")
ax.plot(t, E_ivp_LSODA, label = "solve_ivp: LSODA")
ax.legend()

ax.set_ylim(-0.25, -0)
ax.set_title("Energy vs time (~20,000 periods) for various integrators")
plt.show()