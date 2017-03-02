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
                                              bbox_x[0],
                                              bbox_x[1],
                                              bbox_y[0],
                                              bbox_y[1])

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

    return Q_map, sd_map, Q_r, sd_r, filename


def get_snaps(directory = "."):
    n_snaps = 0

    while "snapshot_{:03d}.hdf5".format(n_snaps) in os.listdir(directory):
        n_snaps += 1
        if n_snaps > 10000:
            break

    return n_snaps



if __name__ == "__main__":
    # Run in script mode

    import sys
    import os
    from multiprocessing import Pool
    from functools import partial

    if "--save" in sys.argv:
        import pickle

    # Physics Setup
    bbox_x = [-100, 100]
    bbox_y = bbox_x
    res_elem = 5.
    res = survis.helper.get_res(res_elem, bbox_x, bbox_y)

    # Computing Setup
    n_cpus = os.cpu_count() if os.cpu_count() <= 8 else 8

    n_snaps = get_snaps()
    filenames = ["snapshot_{:03d}.hdf5".format(x) for x in range(n_snaps)]

    if "--test" in sys.argv:
        filenames = ['test_data.hdf5']

    # We need more arguments than Pool.map allows
    mapped_process = partial(processing_run, res=res, bbox_x=bbox_x, bbox_y=bbox_y)

    with Pool(processes=n_cpus) as processing_pool:
        result = processing_pool.map(mapped_process, filenames)
