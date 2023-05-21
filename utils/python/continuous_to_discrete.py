"""
    This script was used to convert the continuous time state space model to a discrete state space model.

    Linear model parameters were obtained from reference:
        On the coupling of model predictive control and robust Kalman filtering, by Alberto Zenere and Mattia Zorzi
"""

# Importing necessary libraries
import numpy as np
import scipy.linalg


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


# Define continuous time state space model

## Define A
Ac = np.array([
    [0, 1, 0, 0],
    [-(k_theta) / J_l, -beta_l / J_l, (k_theta) / (rho * J_l), 0],
    [0, 0, 0, 1],
    [(k_theta) / (rho * J_m), 0, (-k_theta) / ((rho**2) * J_m),  - ((Kt**2) / (J_m * R)) - (beta_m / J_m)]
])

## Define B
Bc = np.array([
    [0],
    [0],
    [0],
    [Kt / (J_m * R)]
])

## Define C
C = np.array([
    [1, 0, 0, 0]
])


# Convert the model to discrete-time
sampling_time = 0.1
Ad = scipy.linalg.expm(Ac * sampling_time)
Bd = np.matmul(np.linalg.pinv(Ac), np.matmul((scipy.linalg.expm(Ac * sampling_time)-np.eye(4, 4)), Bc))

print("Discrete time state space: ")
print(f"A:\n {Ad}\n")
print(f"B:\n {Bd}\n")
print(f"C:\n {C}")