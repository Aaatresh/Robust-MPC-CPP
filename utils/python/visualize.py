"""
    This script is to post process the input, output and time stamp data, and visualize them on graphs.
    This visualization script is controller independent.

"""

# Import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser("Visualization script for MPC simulation and testing.")
parser.add_argument("--setpoint", default=None, required=True, type=float, help="Set point of tracking problem")
parser.add_argument("--datafile", default=None, required=True, type=str, help="Path and name of data file to be visualized")
args = parser.parse_args()

# Load data from data file
data = np.loadtxt(args.datafile, delimiter=", ", skiprows=1)

# Extract time, input and output
t, u, y = data[:, 0], data[:, 1], data[:, 2]

# Visualization
fig, ax = plt.subplots(2, 1, figsize=(19.2, 10.8), sharex=True)
fig.suptitle("Plots of Input and Output versus Time")

## Input plot
ax[0].plot(t, u, label="Input")
ax[0].set_title("Input voltage vs Time")
ax[0].set_ylabel("Input Voltage (V)")

## Output plot
ax[1].plot(t, y, label="Output")
ax[1].set_title("Output Load Angle vs Time")
ax[1].axhline(y=args.setpoint, linestyle='--', c='k', label=f"Set point = {args.setpoint} rad")
ax[1].set_ylabel("Load Angle (rad)")
ax[1].set_xlabel("Time (s)")

plt.legend()
plt.show()