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


def sound_speed_sne(density, f=0.4, F=0.477, fg=0.1, P=300000, G=4.302e-6):
    """ Returns the sound speed for the supernovae driven model in Martizzi 2015
        which has a constant entropy. """

    entropy = 4.5 * (f/F)**(3./2.) * G**(3./4.) * P**(1./2.) * fg **(-1)

    return np.sqrt((5./4.) * entropy * density**(1./4.))


def Q_gas(sound_speed, kappa, density, surface_density, G=4.302e-6):
    # G given in kpc/msun kms^2
    c_s = sound_speed(density)
    c_s[surface_density == 0] = 0
    sd_masked = surface_density + (surface_density == 0)

    return ((c_s * kappa)/(np.pi * G * sd_masked))


def Q_star():
    return
