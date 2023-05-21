"""
Performance of the robust MPC when the approximated model differs from the actual model
"""

# Importing necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import time

# Linear model parameters
L = 0  # armature coil inductance
# J_m = 0.5  # motor inertia
J_m = 0.5  # motor inertia
beta_m = 0.1   # motor viscous friction coefficient
R = 20   # resistance of armature
Kt = 10  # motor constant
rho = 20  # gear ratio
k_theta = 1280.2  # torsional rigidity
J_l = 25  # nominal load inertia
# J_l = 10  # nominal load inertia
beta_l = 25  # load viscous friction coefficient

del_t = 0.1

# Define A
Ac = np.array([
    [0, 1, 0, 0],
    [-(k_theta) / J_l, -beta_l / J_l, (k_theta) / (rho * J_l), 0],
    [0, 0, 0, 1],
    [(k_theta) / (rho * J_m), 0, (-k_theta) / ((rho**2) * J_m),  - ((Kt**2) / (J_m * R)) - (beta_m / J_m)]
])

# Define B
Bc = np.array([
    [0],
    [0],
    [0],
    [Kt / (J_m * R)]
])

# Define C
C = np.array([
    [1, 0, 0, 0]
])

import scipy.linalg

# Convert the model to discrete-time
T = 0.1
A = scipy.linalg.expm(Ac*T)
B = np.matmul(np.linalg.pinv(Ac), np.matmul((scipy.linalg.expm(Ac*T)-np.eye(4, 4)), Bc))

# prediction and control horizons
Hp = 10
Hu = 3

# Set point
r = np.pi / 2
rt = r * np.ones((Hp, 1))

# Weight matrices Rk and Qk
Rk = 0.01
Rk_diag = Rk * np.ones((Hu,))
Rk_mat = np.diag(Rk_diag)

Qk = 5e3
Qk_diag = Qk * np.ones((Hp,))
Qk_mat = np.diag(Qk_diag)

# This is used to extract the first time step's control
W = np.array([
    [1, 0, 0]
])

# Making psi
psi_T = np.matmul(C, A).T
for c in range(2, Hp+1):
    psi_T = np.hstack((psi_T, np.matmul(C, np.linalg.matrix_power(A, c)).T))

psi = psi_T.T

# Making gamma
gamma = 1 * np.ones((Hp, Hu))
for col in range(Hu):

    col_data = []

    for z in range(col):
        col_data.append(np.array([[0]]))

    for row in range(Hp-col):
        col_data.append(np.matmul(C, np.matmul(np.linalg.matrix_power(A, row), B)))

    gamma[:, col] = np.array(col_data).squeeze()


def get_utt(xtt):
    """
       Solve the unconstrained MPC problem through a closed form expression to find the control input
       given the state estimate at that time step.
    """

    term1 = np.matmul(gamma.T, np.matmul(Qk_mat, gamma)) + Rk_mat

    term2 = np.matmul(gamma.T, np.matmul(Qk_mat, (rt - np.matmul(psi, xtt))))

    utt = np.matmul(W, np.matmul(np.linalg.inv(term1), term2))

    return utt


def eval_kld(Pt, param_t, kld_thresh):
    """
        Evaluate KL divergence given the state's covariance matrix, inverse of the lagrange multiplier and
        the radius of a ball around the nominal model.
    """

    # term1 = -np.log(np.linalg.det(np.linalg.inv(np.eye(2, 2) - (param_t * Pt))))
    term1 = np.log(np.linalg.det(np.eye(4, 4) - (param_t * Pt)))
    term2 = np.trace(np.linalg.inv(np.eye(4, 4) - (param_t * Pt)) - np.eye(4, 4))

    term_total = term1 + term2 - kld_thresh

    return term_total


def bijection_algo(Pt):
    """ Root finding algorithm - Bijection algorithm, given the state's covariance matrix """

    eps = 1e-7
    param1 = eps

    eigs = np.linalg.eigvals(Pt + 1e-8 * np.eye(4, 4))

    lam = np.max(eigs)

    if(((1 / lam) > 2 * eps) and (lam > 1e-5)):
        param2 = (1 / lam) - eps
    else:
        param2 = eps

    param_t = param1

    while(np.abs(param1 - param2) > eps):

        param_t = (param1 + param2) / 2
        kld = eval_kld(Pt, param_t, kld_thresh)

        if(kld < 0):
            param1 = param_t
        else:
            param2 = param_t

    return param_t


def non_lin_dyn(x, u, k):

    """
        Simulate discrete time non-linear dynamics. Find next state given current state, input and time step number.
    """

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

    xdot = np.zeros_like(x)

    xdot[0, 0] = x[1, 0]

    Tfl = (alpha_l0 + (alpha_l1 * np.exp(-alpha_l2 * np.abs(x[1, 0])))) * sign(x[1, 0])
    Ts = (k_theta / rho) * ((x[2, 0] / rho) - x[0, 0])
    xdot[1, 0] = (1 / J_l) * ((rho * Ts) - (beta_l * x[1, 0]) - Tfl)

    xdot[2, 0] = x[3, 0]

    Tfm = (alpha_m0 + (alpha_m1 * np.exp(-alpha_m2 * np.abs(x[3, 0])))) * sign(x[3, 0])
    Im = ((u - (Kt * x[3, 0])) / R) * (1 - np.exp(-R * k * T / L))
    Tm = Kt * Im
    xdot[3, 0] = (1 / J_m) * (Tm - Ts - (beta_m * x[3, 0]) - Tfm)

    return xdot


def sign(x):
    """
        Signum function.
    """
    if(x > 0):
        return 1
    elif(x < 0):
        return -1
    else:
        return 0


# Threshold for KL divergence. In the reference literature, this is mentioned as 'c'
kld_thresh = 0.1

# Initialize yt
y0 = np.array([[0]])
next_yt = y0

# Simulation time parameters
tspan = [0, 20]
samp_time = 0.1

# Initial state covariance and mean
Pt = np.eye(4)
# xtt_1 = np.zeros((4, 1))
xtt_1 = 1e-2 * np.random.randn(4, 1)

all_Ys = []
all_Us = []
all_covs = []
# exit()

# Standard deviation of state and measurement noise
G1 = 1e-2 * np.diag(np.array([0, 1, 1, 1]))
D1 = 1e-2 * np.array([[1, 0, 0, 0]])

t_array = np.arange(tspan[0], tspan[1], samp_time)

# Simulation loop
start_time = time.time()
for e, t in enumerate(np.arange(tspan[0], tspan[1], samp_time)):

    # collect new data
    yt = next_yt

    all_Ys.append(yt[0, 0])
    all_covs.append(Pt[0, 0])

    # Find param_t using the bijection algo
    param_t = bijection_algo(Pt)

    # Determine Vt
    Vt = np.linalg.pinv(np.linalg.pinv(Pt) - (param_t * np.eye(4, 4)))

    # Find Lt
    term1 = np.matmul(Vt, C.T)
    term2 = np.linalg.pinv(np.matmul(C, np.matmul(Vt, C.T)) + np.matmul(D1, D1.T))
    Lt = np.matmul(term1, term2)

    # Prediction step to find xtt
    xtt = xtt_1 + (Lt * (yt - np.matmul(C, xtt_1)))

    # Get control
    utt = get_utt(xtt)

    # temp
    if(utt > 220):
        utt = np.array([[220]])
    elif(utt < -220):
        utt = np.array([[-220]])

    all_Us.append(utt[0, 0])

    # Apply utt to the system
    next_xtt = xtt + (T * non_lin_dyn(xtt, utt, e + 1))
    next_yt = np.matmul(C, next_xtt)

    # Kalman gain
    kg = np.matmul(A, Lt)

    # Update pt
    term1 = np.matmul(A, np.matmul(Vt, A.T))
    term2 = np.matmul(C, np.matmul(Vt, C.T)) + np.matmul(D1, D1.T)
    term3 = np.matmul(kg, np.matmul(term2, kg.T))
    Pt = term1 - term3 + np.matmul(G1, G1.T)

    # Prediction step
    xtt_1 = np.matmul(A, xtt_1) + (kg * (yt - np.matmul(C, xtt_1))) + (B * utt)

end_time = time.time()

print("Execution time: ", (end_time - start_time) / t_array.size)

plt.figure()
plt.plot(t_array, all_Ys, label='Output')
plt.axhline(y=r, color='k', linestyle='--', label='Set-point= (pi / 2)')
plt.fill_between(t_array, y1=[y + np.sqrt(c) for y, c in zip(all_Ys, all_covs)], y2=[y - np.sqrt(c) for y, c in zip(all_Ys, all_covs)], alpha=0.5)
plt.title("Plot of output versus time")
plt.ylabel("Angle (rad)")
plt.xlabel("Time (sec)")
plt.legend()

plt.figure()
plt.plot(all_Us)
plt.title("Plot of control input versus time")
plt.ylabel("Input voltage (V)")
plt.xlabel("Time (sec)")

# Saving all the arrays
# np.save("npy_files/exp2/rob_mpc_adv_dyn_y", np.array(all_Ys))
# np.save("npy_files/exp2/rob_mpc_adv_dyn_u", np.array(all_Us))
# np.save("npy_files/exp2/rob_mpc_adv_dyn_cov", np.array(all_covs))

plt.show()