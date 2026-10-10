import scipy as sp

# Given

P_inf = 20 #kPa
T_inf = 300 #K
Po2 = 260 #kPa
R = 287 #J/(kg*K)
gamma = 1.4
Cp = 1005 #J/(kg*K)

# PART i

from src.normal_shocks import p2_p1
from src.normal_shocks import mach2
from src.isentropic_relations import po_p


# Define function to solve numerically

def f(m1):
    return Po2/P_inf - p2_p1(gamma, m1) * po_p(gamma, mach2(gamma, m1))

# Use Brent's method to solve for M1...

from scipy.optimize import brentq

M1 = brentq(f, 1.2, 4)

print ("M1 =", M1)

# PART ii

# Required stagnation pressure is just the stagnation pressure before the shock (p0_1)

Po_1 = po_p(gamma, M1) * P_inf

print ("Po_1 =", Po_1)

# PART iii

# Required stagnation temperature is just stagnation temperature before the shock

from src.isentropic_relations import stagnation_temperature
from src.isentropic_relations import speed_of_sound

To_1 = stagnation_temperature (T_inf, speed_of_sound(gamma, R, T_inf) * M1, Cp)

print ("To_1 =", To_1)

# PART iv

# This is given in the problem

P_1 = P_inf

print("P_1 =", P_1)

# PART v

# This is the same as To_1 because shocks are adiabatic

To_2 = To_1

print("T_o2 =", To_2)

# PART vi

U_inf = M1 * speed_of_sound(gamma, R, T_inf)

print("U_inf =", U_inf)













