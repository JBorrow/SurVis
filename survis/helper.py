""" Contains helper functions for common plotting operations,
    such as finding the appropriate number of bins given a
    smoothing length. """

import survis.toomre as toom
import survis.fiducial as fid
import numpy as np


def get_res(smoothing, bbox_x, bbox_y):
    """ Smoothing should be a number, bbox_x = [minx, maxx] """

    def res(bbox, sm):
        return int((bbox_x[1] - bbox_x[0])/sm)

    return [res(bbox_x, smoothing), res(bbox_y, smoothing)]


def get_toomre_Q(DG, sound, res_elem):
    """ DG is the datagrid, and sound is a callable function
        that returns the sound speed (see toomre.py).

        res_elem gives the size of the resolution element in simulation units
        (used to convert the mass in a pixel to surface density)

        It automatically masks regions where there are no particles,
        which can be used through cmap.set_bad()."""

    area = res_elem**2
    gas_sd = DG.gas_data['masses']/area
    star_sd = DG.star_data['masses']/area

    gas_v = DG.gas_data['velocities']  # note this is actually v/r
    gas_d = DG.gas_data['densities']

    # Surface density reasoning see 1503.07873v1
    gas_q = toom.Q_gas(sound, gas_v, gas_d, gas_sd + (2./3.)*star_sd)

    # now the sections with no particles must be masked

    return np.ma.array(gas_q, mask=(gas_q == 0.))


def toomre_Q_r(DG, sound, res_elem, max_radius):
    """ Plots the toomre Q as a function of R using fiducial.toomre_Q_gas. """
    # Yes, this is very slow.
    radii = np.arange(res_elem, max_radius, res_elem)
    toomre_Q = []

    for rad in radii:
        toomre_Q.append(fid.toomre_Q_gas(DG, rad, res_elem, sound))

    return toomre_Q


def sd_r(DG, res_elem, max_radius, errors=False):
    """ Plots the surface density as a function of R using
        fiducial.surface_density. """
    # Yes, this is very slow.
    radii = np.arange(res_elem, max_radius, res_elem)
    sd = []

    for rad in radii:
        sd.append(fid.surface_density(DG, rad, res_elem, errors))

    return sd

def n_particles_bins(DG, bins=[0, 0.5, 3, 10, 100]):
    """ Finds the number of particles within the bin radii, useful for seeing
        how the disk stabalises (does it transport mass into the centre?) """
    radii = fid.rss(DG.gas['Coordinates'][()])
    hist, bin_edges = np.histogram(radii, bins)

    return hist, bin_edges

