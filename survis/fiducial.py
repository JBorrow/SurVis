""" Contains extraction functions for things at a fiducial radius over time """

import numpy as np
from survis.toomre import sound_speed

def rms(x):
    return np.sqrt(np.mean(np.square(x), 1))


def rss(x):
    return np.sqrt(np.sum(np.square(x), 1))


def surface_density(DG, R, dR):
    """ Takes the DataGrid, some radius R to find the surface density at over
        some smoothing dR (finds particles within R - dR/2 and R + dR/2) """

    def sd_per_type(data, mass):
        """ Calculates the surface density for a given GADGET dataset and
            particle mass. """

        radii = rss(data['Coordinates'][()])
        radii_mask = np.logical_or((radii < (R - dR)), (radii > (R + dR)))
        within_bounds = np.ma.array(radii, mask=radii_mask)

        n_particles = within_bounds.count()
        m_particles = n_particles * mass

        area_enclosed = 2 * np.pi * R * dR

        return m_particles / area_enclosed

    sd_gas = sd_per_type(DG.gas, DG.gas_mass)
    sd_star = sd_per_type(DG.star, DG.star_mass)

    return [sd_gas, sd_star]


def toomre_Q_gas(DG, R, dR, sound_speed=sound_speed, G=4.302e-6):
    """ Similarly to the above surface_density function, this takes a DataGrid,
        and finds the average toomre Q for the gas within some radius. """

    data = DG.gas

    radii = rss(data['Coordinates'][()])
    radii_mask = np.logical_or((radii < (R - dR)), (radii > (R + dR)))
    vector_mask = np.stack([radii_mask, radii_mask, radii_mask])

    vels = np.mean(rss(np.ma.array(data['Velocities'][()], mask=vector_mask)))
    densities = np.mean(np.ma.array(data['Density'], mask=radii_mask))

    surf_dens = sum(surface_density(DG, R, dR))

    Q = (sound_speed(densities) * (vels/R))/(np.pi * G * surf_dens)

    return Q
