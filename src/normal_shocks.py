def mach2(gamma, m1):
    """Calculate Mach number after a normal shock"""
    return ((1 + 0.5 * (gamma - 1) * m1**2) / (gamma * m1**2 - 0.5 * (gamma - 1)))**0.5

def p2_p1 (gamma, m1):
    """Calculate static pressure ratio P2/P1 across a normal shock"""
    return 1 + ((2 * gamma) / (gamma + 1)) * (m1**2 - 1)

def po2_po1(gamma, M1):
    """Calculate stagnation pressure ratio P02/P02 across a normal shock"""
    return ((((gamma + 1) * 0.5 * M1**2) / (1 + (gamma - 1) * 0.5 * M1**2))**(gamma / (gamma - 1)) *
            (((2 * gamma) / (gamma + 1)) * M1**2 - (gamma - 1) / (gamma +1)) ** (-1 / (gamma - 1)))

def shock_location (Ae, At, Pb_Po1, gamma):
    """Calculate the location of a shock in a choked C-D nozzle"""

    def f(A1):
        from src.quasi_1D_flows import area_mach_relation
        M1 = area_mach_relation(A1, At, gamma, "supersonic")
        Po2_Po1 = po2_po1(gamma, M1)
        At_A2star = Po2_Po1
        Ae_A2star = (Ae / At) * At_A2star
        A2star = Ae / Ae_A2star
        Me = area_mach_relation(Ae, A2star, gamma, "subsonic")
        from src.isentropic_relations import po_p
        Pe_Po2 = 1 / po_p(gamma, Me)
        Pe_Po1 = Pe_Po2 * Po2_Po1
        return Pe_Po1 - Pb_Po1

    from scipy.optimize import brentq
    A1 = brentq(f, At, Ae)

    return A1

test = shock_location(.3198, .1975, 0.784, 1.4)

print(test)

M2 = mach2(1.4, 1.587)
print(M2)