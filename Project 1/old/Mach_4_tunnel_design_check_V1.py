import math as m

# Given (can't change)

gamma = 1.25
L_test = 20 # (cm)
M_e = 4
d_throat_min = 5 # (mm)
P_b_min = 10 # (kPa)

# Variables

P_b = 10 # (kPa)
d_throat = 20 # (mm)
altitude = 21.5 # (km)
alpha = 20 # (deg)

# What is the exit pressure and temperature??

from src.standard_atmosphere import T_std_atm
from src.standard_atmosphere import P_std_atm

T_e = T_std_atm(altitude * 1000)
P_e = P_std_atm(altitude * 1000, T_e)

# Area ratio given by:

from src.isentropic_relations import F

A_Astar = F(gamma, M_e)

# What is the stagnation pressure and temperature?

from src.isentropic_relations import po_p
from src.isentropic_relations import to_t

P_o = po_p(gamma, M_e) * P_e
T_o = to_t(gamma, M_e) * T_e

# Check exit condition

if P_e < P_b:
    print("Overexpanded")
    print("Pe =", P_e, "Pb =", P_b)
    print("P_b/P_e =", P_b/P_e)

# Need to check if a shock can produce this pressure ratio. Limiting case: normal shock

from src.normal_shocks import p2_p1

pb_pe_max = p2_p1(gamma, M_e)

if pb_pe_max < (P_b / P_e):
    print("Pressure ratio is too high. Decrease back pressure")
    print("Pb/Pe_max = ", pb_pe_max)
else:
    print("Shock solution exists for design pressure ratio")
    print("Pb/Pe_max = ", pb_pe_max)
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

d_e = m.sqrt(A_Astar * d_throat**2) # (outputs d_e in mm)

# Check required angle phi to satisfy L_test

phi = m.atan(((d_e / 10) / (2 * L_test))) # (d_e and L_test in cm)

# Required beta is calculated from theta and phi

beta_req = phi + alpha/2

if beta <= beta_req:
    print("Test length requirement satisfied")
    print("beta =", m.degrees(beta), "beta_req =", m.degrees(beta_req))
elif beta > beta_req:
    print("Shock angle is too large")
    print("beta =", m.degrees(beta), "beta_req =", m.degrees(beta_req))

# Check if turning angle satisfies design constraint

alpha = 2 * theta

if m.degrees(alpha) < alpha_max:
    print("Divergence angle design requirement met")
    print("theta =", m.degrees(theta))
else:
    print("Divergence angle is too large. If minimum back pressure is met, reduce the target altitude")
    print("theta =", m.degrees(theta))

