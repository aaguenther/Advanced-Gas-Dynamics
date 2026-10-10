import matplotlib.pyplot as plt
import scipy as scp
import numpy as np

# Given

R = 300 # (J/(kg*K))
gamma = 1.4
U_inf = 600 # (m/s)
T_inf = 100 # (K)
cp = (R * gamma) / (gamma - 1)

# This is a range of the non-dimensional variable eta = y/delta to plot over

eta = np.arange(0, 1.01, 0.01)

# Velocity in boundary layer (from problem)

def u(eta):
    return U_inf * np.sin(eta * (np.pi / 2))

# Adiabatic wall temperature (from problem)

T_aw = T_inf + U_inf**2 / (2 * cp)

# Temperature in boundary layer (from problem)

def bl_static_temp(t_wall, eta):
     return t_wall * (1 + ((T_aw - t_wall) / t_wall) * (u(eta) / U_inf) - u(eta)**2 / (2 * cp * t_wall))

# PART iv

T_wall = np.array([25, 285, T_aw])

# Import functions I made

from src.normal_shocks import mach2
from src.normal_shocks import p2_p1
from src.isentropic_relations import po_p

# Define function for pressure ratio in the boundary layer

def bl_stagnation_pressure(eta, t_wall):
    """Includes logic for calculating pressure ratio based on subsonic/supersonic conditions"""
    mach = u(eta) / (gamma * R * bl_static_temp(t_wall, eta))**0.5
    print(mach)
    if mach >= 1:
        return p2_p1(gamma, mach) * po_p(gamma, mach2(gamma, mach))
    else:
        return po_p(gamma, mach)

# Create an array with pressure ratios in the boundary layer at each value of T_wall

Po = np.array([
    [bl_stagnation_pressure(e, t) for e in eta]
    for t in T_wall
])

# Create the plot

plt.plot(Po[0], eta)
plt.plot(Po[1], eta)
plt.plot(Po[2], eta)

# Format plot

plt.ticklabel_format(axis='x', style='plain', useOffset=False)
plt.xlabel("Po/P_inf")
plt.ylabel("y/delta")
plt.title("Stagnation Pressure Normalized w/ Freestream Pressure")
plt.legend(["25 K", "285 K", "T_aw"])

plt.grid()

# Save the plot

plt.savefig("plots/Po_P_inf.png")
plt.show()

# PART v

# This is the same as part iv since P = P_inf in a boundary layer

# Create the plot

plt.plot(Po[0], eta)
plt.plot(Po[1], eta)
plt.plot(Po[2], eta)

# Format plot

plt.ticklabel_format(axis='x', style='plain', useOffset=False)
plt.xlabel("Po/P")
plt.ylabel("y/delta")
plt.title("Stagnation Pressure Normalized w/ Static Pressure")
plt.legend(["25 K", "285 K", "T_aw"])

plt.grid()

# Save the plot

plt.savefig("plots/Po_P.png")
plt.show()

# PART vi

# Mach number in boundary layer

def bl_mach(eta, t_wall):
    return u(eta) / (gamma * R * bl_static_temp(t_wall, eta))**0.5

M = np.array([
    [bl_mach(e, t) for e in eta]
    for t in T_wall
    ])

# Create the plot

plt.plot(M[0], eta)
plt.plot(M[1], eta)
plt.plot(M[2], eta)

# Format plot

plt.ticklabel_format(axis='x', style='plain', useOffset=False)
plt.xlabel("M")
plt.ylabel("y/delta")
plt.title("Mach Number in Supersonic Boundary Layer")
plt.legend(["25 K", "285 K", "T_aw"])

plt.grid()

# Save the plot

plt.savefig("plots/mach_number.png")
plt.show()

# PART vii

T_T_inf = np.array([
    [bl_static_temp(t, e) / T_inf for e in eta]
    for t in T_wall
    ])

# Create the plot

plt.plot(T_T_inf[0], eta)
plt.plot(T_T_inf[1], eta)
plt.plot(T_T_inf[2], eta)

# Format plot

plt.ticklabel_format(axis='x', style='plain', useOffset=False)
plt.xlabel("T/T_inf")
plt.ylabel("y/delta")
plt.title("Temperature Normalized by Freestream Temperature")
plt.legend(["25 K", "285 K", "T_aw"])

plt.grid()

# Save the plot

plt.savefig("plots/Temperature.png")
plt.show()

