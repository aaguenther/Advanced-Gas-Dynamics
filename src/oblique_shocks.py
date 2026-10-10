import numpy as np
from scipy.optimize import brentq
from scipy.optimize import minimize_scalar

# Run time on this one is pretty bad...

def theta_beta_m(m1, theta_deg, gamma):
    """theta_beta_m relation for oblique shocks that returns shock angle (beta)"""
    theta = np.radians(theta_deg)
    # Residual function to iterate over
    def f(beta):
        return (
            np.tan(theta)
            - 2 * (1 / np.tan(beta))
            * ((m1**2 * np.sin(beta)**2 - 1)
               / (m1**2 * (gamma + np.cos(2 * beta)) + 2))
        )

    # Find mach angle
    mu = np.arcsin(1 / m1)

    # Find beta for the maximum turning angle for weak solution using minimize_scalar (easier than using first derivative of theta(beta))
    result = minimize_scalar(
        # The maximum of theta(beta) is equal to the minimum of -theta(beta)
        lambda beta: -(
            np.arctan(
                2 / np.tan(beta)
                * (m1 ** 2 * np.sin(beta) ** 2 - 1)
                / (m1 ** 2 * (gamma + np.cos(2 * beta)) + 2)
            )
        ),
        bounds= (mu, np.pi / 2),
        method='bounded'
    )

    beta_max = result.x

    # Weak-shock solution
    beta = brentq(f, mu, beta_max)

    return np.rad2deg(beta)

def m1_to_m2(m1, theta, beta, gamma):
    """Calculates mach number change across an oblique shock"""
    m1n = np.sin(beta) * m1
    m2n = np.sqrt((1 + 0.5 * (gamma - 1) * m1n**2) /
                  (gamma * m1n**2 - 0.5 * (gamma - 1)))
    m2 = m2n / np.sin(beta - theta)
    return m2

def p2_p1_oblique(m1, beta, gamma):
    """Calculates pressure ratio across an oblique shock"""
    m1n = np.sin(np.deg2rad(beta)) * m1
    from src.normal_shocks import p2_p1
    return p2_p1(gamma, m1n)


