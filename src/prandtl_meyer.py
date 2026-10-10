import numpy as np

def prandtl_meyer(M, gamma):
    nu = (
        np.sqrt((gamma + 1) / (gamma - 1))
        * np.arctan(
            np.sqrt(
                ((gamma - 1) / (gamma + 1))
                * (M**2 - 1)
            )
        )
        - np.arctan(np.sqrt(M**2 - 1))
    )

    return np.rad2deg(nu)

