""" Does some common tasks, including

    + Creating a movie of toomre Q against time (snapshot number)
    + Creating a movie of surface density against time (snapshot number)
    + Creating a plot of surface density, Q at the solar radius against time
    + Creating a movie of the surface density evolution map
    + Creating a movie of the toomre Q evolution map

    By using the argument --save, one can also save these things to a file
    called processed_variables.pkl

    Please note that by default this uses the Supernovae equation of state
    rather than an isothermal one when calculting the Toomre Q parameter.
"""

import numpy as np
import matplotlib.pyplot as plt
import survis

# Constants
solar_radius = 8.  # kpc
smoothing = 0.2 * 2  # kpc

def processing_run(filename, res, bbox_x, bbox_y):
    """ Generates the processed data out of the snapshot """

    data_grid = survis.preprocess.DataGridder(filename,
                                              res[0],
                                              res[1],
                                              bbox_x,
                                              bbox_y)

    # First we extract the maps
    Q_map = survis.helper.get_toomre_Q(data_grid,
                                       survis.toomre.sound_speed_sne,
                                       res[0])

    # Normally the masses of each element are given, we must divide by size
    sd_map = data_grid.gas_data['masses']/(res[0]*res[1])


    # Now the values at a given radius
    sd_r = survis.fiducial.surface_density(data_grid, solar_radius, smoothing)
    Q_r = survis.fiducial.toomre_Q_gas(data_grid,
                                       solar_radius,
                                       smoothing,
                                       survis.toomre.sound_speed_sne)

    return Q_map, sd_map, Q_r, sd_r



if __name__ == "__main__"
    # Run in script mode

    import sys
    import os

    from multiprocessing import Pool

    while "snapshot_{:03d}.hdf5".format(n_snaps) in os.listdir("."):
        n_snaps += 1
        if n_snaps > 10000:
            break
