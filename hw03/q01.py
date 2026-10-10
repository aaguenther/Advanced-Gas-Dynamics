import math
import numpy as np
import matplotlib.pyplot as plt

# GIVEN

A_hole = 0.01 # m^2
A_airlock = math.pi # m^2
P_3 = 0 # Pa
P_1 = 100E3 # Pa
T_1 = 300 # K
gamma = 5/3
R = 300 # J/kg*K

# PART i

from src.quasi_1D_flows import area_mach_relation

M_airlock = area_mach_relation(A_airlock, A_hole, gamma, 'subsonic')

print("Mach number upstream in airlock =", M_airlock)

# PART ii

from src.isentropic_relations import po_p
from src.isentropic_relations import to_t
from src.isentropic_relations import speed_of_sound

P_airlock = (1 / po_p(gamma, M_airlock)) * P_1
T_airlock = (1 / to_t(gamma, M_airlock)) * T_1
U_airlock = M_airlock * speed_of_sound(gamma, R, T_airlock)

print("Pressure upstream in airlock (Pa) =", P_airlock)
print("Temperature upstream in airlock (K) =", T_airlock)
print("Velocity upstream in airlock (m/s) =", U_airlock)

# PART iv

# Array of specific heat ratios

gamma_range = np.arange(1.05, 5/3, 0.01)

# Initialize array of mach numbers

mach_numbers = []

# Iterate over gammas to fill array

for i in gamma_range:
    M_airlock = area_mach_relation(A_airlock, A_hole, i, 'subsonic')
    mach_numbers.append(M_airlock)

# Plot Mach numbers vs. gamma

plt.plot(gamma_range, mach_numbers)
plt.xlabel("gamma")
plt.ylabel("M_airlock")
plt.title("Mach Number in Airlock vs. Gamma")
plt.show()

# PART v

from src.isentropic_relations import rho_o_rho

# Calculate total density using ideal gas EoS

rho_1 = P_1 / (R * T_1)

print('Total density (kg/m^3) =', rho_1)

# Calculate sonic conditions at the hole

rho_star = (1 / rho_o_rho(gamma, 1)) * rho_1
T_star = (1 / to_t(gamma, 1)) * T_1
P_star = (1 / po_p(gamma, 1)) * P_1

print('Density at hole (kg/m^3) =', rho_star)
print('Temperature at hole (K) =', T_star)
print('Pressure at hole (Pa) =', P_star)

# Calculate choked mass flow rate

m_dot_choked = rho_star * speed_of_sound(gamma, R, T_star) * A_hole

print('Choked mass flow rate (kg/s) =', m_dot_choked)




