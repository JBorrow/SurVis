# Contains routines that visualise the toomre Q parameter.
import numpy as np

def sound_speed(density):
    # Normally contains information about the equation of state, as this is
    # formally dP/d\rho.
    gamma = 5./3.
    R = 8.314
    M = 0.001
    T = 1e4
    original = np.ones_like(density) * np.sqrt((gamma * R * T)/(M)) # m/s
    return original * 0.001 # km/s


def Q_gas(sound_speed, kappa, density, surface_density, G=4.302e-6):
    # G given in kpc/msun kms^2
    c_s = sound_speed(density)
    c_s[surface_density == 0] = 0
    sd_masked = surface_density + (surface_density == 0)

    # our surface density is given in msun/kpc^2 so we need a conversion factor
    return 3.086e16 * ((c_s * kappa)/(np.pi * G * sd_masked))


def Q_star():
    return
