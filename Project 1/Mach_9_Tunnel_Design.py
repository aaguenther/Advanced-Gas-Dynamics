import math as m
import numpy as np
import matplotlib.pyplot as plt

# This wind tunnel has the following design parameters:

gamma = 1.25
R = 287 # (J/kg*K)
M_e = 9
P_b = 101.325 # (kPa)
d_throat = 5 # (mm)
altitude = 25 # (km)
T_o_initial = 295 # (K)
alpha = 20 # (deg)
# Geometry in converging section
alpha_c = 30 # degrees
Ac_Astar = 4 # contraction ratio

# Exit pressure and temperature

from src.standard_atmosphere import T_std_atm
from src.standard_atmosphere import P_std_atm

T_e = T_std_atm(altitude * 1000) # Outputs in deg C
P_e = P_std_atm(altitude * 1000, T_e)

# Area ratio given by:

from src.isentropic_relations import F

A_Astar_Exit = F(gamma, M_e)

# What is the stagnation pressure and temperature?

from src.isentropic_relations import po_p
from src.isentropic_relations import to_t

P_o = po_p(gamma, M_e) * P_e
T_o = to_t(gamma, M_e) * (T_e + 273.15) # Outputs K

# Calculate beta for P_b/P_e

beta = m.asin(m.sqrt(
        ((gamma + 1) * (P_b / P_e) + (gamma - 1))
        / (2 * gamma * M_e ** 2)
        ))

# Calculate the flow turning angle

theta = m.atan(
        2 * (1 / m.tan(beta)) *
        ((M_e ** 2 * m.sin(beta) ** 2 - 1) /
         (M_e ** 2 * (gamma + m.cos(2 * beta)) + 2))
         )

# Exit diameter given by:

d_e = m.sqrt(A_Astar_Exit * d_throat**2) # (outputs d_e in mm)
A_e = (m.pi * d_e ** 2) / 4

# Calculate angle phi

phi = beta - m.radians(alpha/2)

# Calculate test section length

L_test = d_e / (2 * m.tan(phi))

# Create plots

# First, need length of nozzle

# Length of diverging section

r_e = d_e / 2
r_throat = d_throat / 2

L_d = (r_e - r_throat) / m.tan(m.radians(alpha) / 2)

# Length of converging section

r_i = m.sqrt(Ac_Astar * r_throat**2)
L_c = (r_i - r_throat) / m.tan(m.radians(alpha_c / 2))

# Total nozzle length

L = L_c + L_d

# Throat area

A_throat = m.pi * r_throat**2 # (mm^2)

# Mass flow rate

A_throat_m2 = A_throat / 1e6
P_o_Pa = P_o * 1000

m_dot = (
    A_throat_m2
    * P_o_Pa
    * m.sqrt(gamma / (R * T_o))
    * (2 / (gamma + 1)) **
      ((gamma + 1) / (2 * (gamma - 1)))
)

# Required heat input

cp = (gamma * R) / (gamma - 1) # (J/kg*K)
Q_dot_in = m_dot * cp * (T_o - T_o_initial) # (W)

# Create array to evaluate nozzle values over

L_array = np.linspace(0, L, 500)

# Initialize arrays for M, A, A/Astar, Po, and To (plotted on whole nozzle)

M = []
A = []
A_Astar = []
P_o_array = []
T_o_array = []

from scipy.optimize import brentq

# Iterate to fill array

for i in range(0, len(L_array)):

    # Converging section

    if L_array[i] <= L_c:
        # Local area and radius
        r = ((r_throat - r_i) / L_c) * L_array[i] + r_i
        A_i = m.pi * r**2
        A_Astar_i = A_i / A_throat
        # Solve F(M) for M using area ratio
        def f(M):
            return A_Astar_i - F(gamma, M)
        M_i = brentq(f, 1e-6, 1 - 1e-6)
        P_o_i = P_o
        T_o_i = T_o

    # Diverging section

    else:
        # Local area and radius
        r = r_throat + ((r_e - r_throat) / L_d) * (L_array[i] - L_c)
        A_i = m.pi * r ** 2
        A_Astar_i = A_i / A_throat
        # Solve F(M) for M using area ratio
        def f(M):
            return A_Astar_i - F(gamma, M)
        M_i = brentq(f, 1 + 1e-6, M_e + 3)
        P_o_i = P_o
        T_o_i = T_o

    # Append values to arrays

    M.append(M_i)
    A.append(A_i)
    A_Astar.append(A_Astar_i)
    P_o_array.append(P_o_i)
    T_o_array.append(T_o_i)

# Initialize arrays for P, T (plotted along diverging section)

P = []
T = []

# Create an array for location and initialize an array for Mach number along diverging section

L_diverging_array = np.linspace(L_c, L_c + L_d, 500)
M_diverging = []

for i in range(0, len(L_diverging_array)):
    # Local area and radius
    r = r_throat + ((r_e - r_throat) / L_d) * (L_diverging_array[i] - L_c)
    A_i = m.pi * r ** 2
    A_Astar_i = A_i / A_throat
    # Solve F(M) for M using area ratio
    def f(M):
        return A_Astar_i - F(gamma, M)
    M_i = brentq(f, 1, M_e + 3)
    P_i = P_o * (1 / po_p(gamma, M_i))
    T_i = T_o * (1/ to_t(gamma, M_i))

    # Append values to arrays

    P.append(P_i)
    T.append(T_i)
    M_diverging.append(M_i)

print("------------------------------------------------------------")
print("---------------------SUMMARY OF RESULTS---------------------")
print("------------------------------------------------------------")
print("Test Conditions")
print("------------------------------------------------------------")
print("Altitude (h) =", altitude, "m")
print("Flight Mach Number (M_e) =", M_e)
print("------------------------------------------------------------")
print("Performance Parameters")
print("------------------------------------------------------------")
print("Test Section Length (L_test) = ", L_test / 10, "cm")
print("Mass Flow Rate (mdot) = ", m_dot, "kg/s")
print("Heat Input Required (Q_dot_in) = ", Q_dot_in / 1000, "kW")
print("------------------------------------------------------------")
print("Nozzle Geometry")
print("------------------------------------------------------------")
print("alpha = ", alpha, "degrees")
print("Throat Diameter (d_throat) = ", d_throat, "mm")
print("Throat Area (A_throat) =", A_throat, "mm^2")
print("Exit Diameter (d_e) = ", d_e, "mm")
print("Exit Area (A_e) =", A_e, "mm^2")
print("Inlet Diameter (d_i) = ", r_i * 2, "mm")
print("Divergence Area Ratio (A_e/A*) =", A_e/A_throat)
print("Diverging section length (L_d) =", L_d, "mm")
print("Converging section length (L_c) =", L_c, "mm")
print("Total Nozzle Length (L) =", L, "mm")
print("------------------------------------------------------------")
print("Shock Geometry")
print("------------------------------------------------------------")
print("Exit flow angle (slip line) (theta) = ", m.degrees(theta), "degrees")
print("Shock angle (beta) = ", m.degrees(beta), "degrees")
print("------------------------------------------------------------")
print("Flow Parameters")
print("------------------------------------------------------------")
print("Stagnation Pressure (P_o) =", P_o, "kPa")
print("Stagnation Temperature (T_o) =", T_o, "K")
print("Exit Pressure (P_e) =", P_e, "kPa")
print("Exit Temperature (T_e) =", T_e + 273.15, "K")
print("Back Pressure (P_b) =", P_b, "kPa")
print("Pressure Ratio Across Oblique Shock (P_b/P_e) =", P_b / P_e)
print("------------------------------------------------------------")

# Create Plots

plt.plot(L_array, M)
plt.title("Mach number vs. Nozzle Axial Location (Mach 9)")
plt.xlabel("Nozzle Axial Location (mm)")
plt.ylabel("Mach Number")
plt.grid()
plt.show()

plt.plot(L_array, A)
plt.title("Area vs. Nozzle Axial Location (Mach 9")
plt.xlabel("Nozzle Axial Location (mm)")
plt.ylabel("Area (mm^2)")
plt.grid()
plt.show()

plt.plot(L_array, A_Astar)
plt.title("A/A* vs. Nozzle Axial Location (Mach 9)")
plt.xlabel("Nozzle Axial Location (mm)")
plt.ylabel("A/A*")
plt.grid()
plt.show()

plt.plot(L_array, P_o_array)
plt.title("Stagnation Pressure vs. Nozzle Axial Location (Mach 9)")
plt.xlabel("Nozzle Axial Location (mm)")
plt.ylabel("Po (kPa)")
plt.grid()
plt.show()

plt.plot(L_array, T_o_array)
plt.title("Stagnation Temperature vs. Nozzle Axial Location (Mach 9)")
plt.xlabel("Nozzle Axial Location (mm)")
plt.ylabel("To (K)")
plt.grid()
plt.show()

plt.plot(L_diverging_array, P)
plt.title("Pressure vs. Nozzle Axial Location Along Diverging Section (Mach 9)")
plt.xlabel("Nozzle Axial Location (mm)")
plt.ylabel("P (kPa)")
plt.grid()
plt.show()

plt.plot(L_diverging_array, T)
plt.title("Temperature vs. Nozzle Axial Location Along Diverging Section (Mach 9)")
plt.xlabel("Nozzle Axial Location (mm)")
plt.ylabel("T (K)")
plt.grid()
plt.show()