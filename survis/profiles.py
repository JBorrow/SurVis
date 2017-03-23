""" This section of the code fits profiles to the disk and returns scale
    lenghts, heights, radii... Whatever is relevant.

    It is useful to compare these to the local Jeans' length to check if
    the assumption of Schaye 2001 is valid. """

import survis.fiducial as fid
import numpy as np

from scipy.optimize import curve_fit


def bin_cent(bins):
    """ We get bin edges from our histogramming, to plot we want
        the centers of the bins! """
    return 0.5*(bins[1:] + bins[:-1])


def radial_profile(DG, bin_width=0.4):
    """ Takes the data grid for a galaxy and fits the profile radially.
        Expects an exponential surface density profile."""
    radii = fid.rss(DG.gas['Coordinates'][()])

    n, bins = np.histogram(radii, bin_width)
    bincenters = bin_cent(bins)

    def to_fit(r, norm, R):
        return r*norm*np.exp(-r/R)

    popt, pcov = curve_fit(to_fit, bincenters, n)

    return popt[1], np.sqrt(pcov[1,1])


def vertical_profile(DG, bin_width=0.2, min=-10, max=10):
    z = DG.gas['Coordinates'][:, 2]
    
    bins = np.arange(min, max, bin_width)
    
    n, bins = np.histogram(z, bins)
    bincenters = bin_cent(bins)

    def to_fit(z, norm, Z):
        return norm*(1/(np.cosh(z/Z)**2))

    popt, pcov = curve_fit(to_fit, bincenters, n)

    return popt[1], np.sqrt(pcov[1,1])


def local_jeans_length(DG):
    return 0
