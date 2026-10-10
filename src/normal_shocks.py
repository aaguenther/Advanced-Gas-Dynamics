def mach2(gamma, m1):
    """Calculate Mach number after a normal shock"""
    return ((1 + 0.5 * (gamma - 1) * m1**2) / (gamma * m1**2 - 0.5 * (gamma - 1)))**0.5

def p2_p1 (gamma, m1):
    """Calculate static pressure ratio P2/P1 across a normal shock"""
    return 1 + ((2 * gamma) / (gamma + 1)) * (m1**2 - 1)



