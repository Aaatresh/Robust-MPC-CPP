"""
Simulating non-linear dynamics
"""

# Importing necessary libraries
import numpy as np
import matplotlib.pyplot as plt

# Model parameters
eps_max = 0
eps_min = 0

L = 0.8  # armature coil inductance
J_m = 0.5 * (1 + eps_max)  # motor inertia
beta_m = 0.1 * (1 + eps_max)  # motor viscous friction coefficient
R = 20 * (1 + eps_min)  # resistance of armature
Kt = 10 * (1 + eps_max)  # motor constant
rho = 20 * (1 + eps_min)  # gear ratio
k_theta = 1280.2 * (1 + eps_min)  # torsional rigidity
J_l = 25 * (1 - eps_max)  # nominal load inertia
beta_l = 25 * (1 + eps_max)  # load viscous friction coefficient
alpha_l0, alpha_l1, alpha_l2 = [0.5, 10, 0.5]  # load non-linear friction parameters
alpha_m0, alpha_m1, alpha_m2 = [0.1, 2, 0.5]  # motor non-linear friction params

# Sampling time
T = 0.01

def sign(x):
    """ The sign function """

    if(x > 0):
        return 1
    elif(x < 0):
        return -1
    else:
        return 0

def f(x, u, k):

    """ Non-linear state transition function """

    xdot = np.zeros_like(x)

    xdot[0, 0] = x[1, 0]

    Tfl = (alpha_l0 + (alpha_l1 * np.exp(-alpha_l2 * np.abs(x[1, 0])))) * sign(x[1, 0])
    Ts = (k_theta / rho) * ((x[2, 0] / rho) - x[0, 0])
    xdot[1, 0] = (1 / J_l) * ((rho * Ts) - (beta_l * x[1, 0]) - Tfl)

    xdot[2, 0] = x[3, 0]

    Tfm = (alpha_m0 + (alpha_m1 * np.exp(-alpha_m2 * np.abs(x[3, 0])))) * sign(x[3, 0])
    Im = ((u - (Kt * x[3, 0])) / R) * (1 - np.exp(-R * k * T / L))
    Tm = Kt * Im
    xdot[3, 0] = (1 / J_m) * (Tm - Ts - (beta_m * x[3 ,0]) - Tfm)

    return xdot

# Define C
C = np.array([
    [1, 0, 0, 0]
])

""" Lets simulate and visualize these dynamics """

x0 = np.zeros((4, 1))
t_array = np.arange(0, 5, T)
# u_array = 100 * np.exp(-np.power(t_array - 2.5, 2))  # Gaussian input
u_array = 100 * np.exp(-2 * t_array)  # exponentially decreasing input
# u_array = np.concatenate((np.ones((t_array.size // 2,)), np.zeros((t_array.size - (t_array.size // 2),))))  # step input
# u_array = np.ones((t_array.size,)) # all ones input


# Visualizing the input
# plt.figure()
# plt.title("Voltage input")
# plt.plot(t_array, u_array)
# plt.show()
# exit()

# Forward propagating the dynamics
xk_1 = x0
Y = []
for e, t_step in enumerate(t_array):
    xk = xk_1 + (T * f(xk_1, u_array[e], e + 1))

    xk_1 = xk

    Y.append(xk[0, 0])

# Visualizing the output and the input
fig, ax = plt.subplots(2, 1)
ax[0].set_title("Load angle")
ax[0].plot(Y)

ax[1].set_title("Voltage input")
ax[1].plot(u_array)

plt.show()
