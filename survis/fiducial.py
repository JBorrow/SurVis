""" Contains extraction functions for things at a fiducial radius over time """

import numpy as np

def surface_density(DG, R, dR):
    """ Takes the DataGrid, some radius R to find the surface density at over
        some smoothing dR (finds particles within R - dR/2 and R + dR/2) """

    def sd_per_type(data, mass):
        """ Calculates the surface density for a given GADGET dataset and
            particle mass. """

        radii = np.sqrt(np.sum(np.square(data['Coordinates'][()]), 1)) - R
        radii_mask = np.logical_or((radii < -dR), (radii > dR))
        within_bounds = np.ma.array(radii, mask=radii_mask)

        n_particles = within_bounds.count()
        m_particles = n_particles * mass

        area_enclosed = 4 * np.pi * R * dR

        return m_particles / area_enclosed

    sd_gas = sd_per_type(DG.gas, DG.gas_mass)
    sd_star = sd_per_type(DG.star, DG.star_mass)

    return [sd_gas, sd_star]


def toomre_Q():
    return 0.
