# Contains routines that visualise the toomre Q parameter.
import numpy as np

def sound_speed(density):
    # Normally contains information about the equation of state, as this is
    # formally dP/d\rho.
    gamma = 5./3.
    R = 8.314
    M = 0.001
    T = 1e4
    return np.ones_like(density) * np.sqrt((gamma * R * T)/(M)) # m/s


def Q_gas(sound_speed, kappa, density, surface_density, G):
    c_s = sound_speed(density)

    return (c_s * kappa)/(np.pi * G * surface_density)
