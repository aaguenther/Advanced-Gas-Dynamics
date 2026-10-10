# FROM https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html

import math as m

def T_std_atm (h):
    """Returns standard atmospheric temperature (Celsius) for a given altitude (meters)"""
    if h < 0:
        print("Not a valid input. h must be greater than 0")
        return
    if h <= 11000:
        T = 15.04 - .00649 * h
    elif 11000 < h < 25000:
        T = -56.46
    elif h >= 25000:
        T = -131.21 + .00299 * h
    return T

def P_std_atm (h, T):
    """Returns standard atmospheric pressure (kPa) for a given altitude (meters)"""
    if h < 0:
        print("Not a valid input. h must be greater than 0")
        return
    if T < -273.15:
        print("Not a valid input. T must be greater than -273.15")
        return
    if h <= 11000:
        P = 101.29 * ((T + 273.1) / 288.08) ** (5.256)
    elif 11000 < h < 25000:
        P = 22.65 * m.exp(1.73 - .000157 * h)
    elif h >= 25000:
        P = 2.488 * ((T + 273.1) / 216.6) ** (-11.388)
    return P


