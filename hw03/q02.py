import numpy as np

# GIVEN

m1 = 5
theta_2 = 15 # deg
theta_3 = 15 # deg
gamma = 1.4

# Solve for region 2 using oblique shock relations

from src.oblique_shocks import theta_beta_m
from src.oblique_shocks import m1_to_m2
from src.oblique_shocks import p2_p1_oblique

beta_2 = theta_beta_m(m1, theta_2, gamma)
m2 = m1_to_m2(m1, theta_2, beta_2, gamma)
p2_p1 = p2_p1_oblique(m1, beta_2, gamma)

print(m2)

# Solve for region 3 using oblique shock relations

beta_3 = theta_beta_m(m2, theta_3, gamma)
m3 = m1_to_m2(m2, theta_3, beta_3, gamma)
p3_p2 = p2_p1_oblique(m2, beta_3, gamma)

print(m3)

# Need pressure ratio wrt incoming flow

p3_p1 = p3_p2 * p2_p1

# Now things get hard.

# This solution involves guessing a value for the slip line angle (call it alpha), then iterating until the pressures
# across the slip line are equal

# If interaction is type 1, the flow region near the wall expands, and the flow across the slip line compresses (ex. 4 and 5)
# If interaction is type 2, the flow region near the wall compresses, and the flow across the slip line expands (ex. 7 and 8)

from src.prandtl_meyer import prandtl_meyer
from scipy.optimize import brentq
from src.isentropic_relations import p2_p1

def slip_line_interaction(m1, m2, p1_pref, p2_pref, gamma, type):
    def f(alpha):
        if type == 1:
            # Oblique shock interaction across the slip line
            beta_3 = theta_beta_m(m1, alpha, gamma)
            p3_pref = p1_pref * p2_p1_oblique(m1, beta_3, gamma)
            # PM expansion near the wall
            delta = alpha - 30
            nu_2 = prandtl_meyer(m2, gamma)
            def f_m4(m4):
                return prandtl_meyer(m4, gamma) - (nu_2 + delta)
            m4 = brentq(f_m4, m2, 5)
            p4_p2 = p2_p1(m4, m2, gamma)
            p4_pref = p4_p2 * p2_pref
            return p4_pref - p3_pref
        else:
            return ("Not a valid input")
    alpha = brentq(f, 30.0001, 45)
    return alpha

alpha_CE = slip_line_interaction(m1, m3, 1, p3_p1, gamma, 1)

print(alpha_CE)






