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
import matplotlib.cm as cm
import matplotlib.colors as col
import matplotlib.animation as animation
from tqdm import tqdm
import survis

# Constants
solar_radius = 8.  # kpc
smoothing = 0.2 * 2  # kpc

def processing_run(filename, res, bbox_x, bbox_y, callback=None):
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

    if not (callback is None):
        callback()

    return Q_map, sd_map, Q_r, sd_r, filename[9:12]


def get_snaps(directory = "."):
    n_snaps = 0

    while "snapshot_{:03d}.hdf5".format(n_snaps) in os.listdir(directory):
        n_snaps += 1
        if n_snaps > 10000:
            break

    return n_snaps


def make_movie_imshow(data, bad_color='black', log=False):
    images = []

    fig = plt.figure()
    colormap = cm.get_cmap('viridis')
    c_scale = col.Normalize(vmin=0 , vmax=3)
    colormap.set_bad(bad_color, -1)


    for item in tqdm(data, desc="Data plotting"):
        images.append([plt.imshow(item, cmap=colormap, vmin=0, vmax=3)])

    fig.colorbar(cmap=colormap, norm=c_scale, mappable=images[-1][0])

    return animation.ArtistAnimation(fig, images, interval=50, repeat_delay=3000,
                                   blit=True)

def make_linear_plot(data, ylabel, ymin=0, ymax=5.):
    n_snaps = len(data)
    fig, ax = plt.subplots()

    ax.plot(np.arange(n_snaps), np.array(data))

    ax.set_xlim(0, n_snaps)
    ax.set_ylim(ymin, ymax)

    ax.set_xlabel("Snapshot number")
    ax.set_ylabel(ylabel)

    return fig, ax

def make_plots(result, make_movies=True):
    result = np.array(result)
    Q_maps = result.T[0]
    sd_maps = result.T[1]
    Q_r = result.T[2]
    # We have to manually convert the lists here from when they get pickled
    sd_r = np.array([np.array(x) for x in result.T[3]]).T
    sd_r_gas = sd_r[0, :]
    sd_r_star = sd_r[1, :]

    plt.imshow(Q_maps[250]/1e5)
    plt.colorbar()
    plt.show()


    Q_fig, Q_ax = make_linear_plot(Q_r, "Toomre $Q$", 0.5, 3.0)
    Q_fig.savefig("Q_fig.pdf")

    print("Writing plots")
    sd_fig, sd_ax = make_linear_plot(sd_r_gas, "Surface Density ($M_\odot$ pc$^{-2}$)", 0, 1e5)
    sd_fig.savefig("sd_gas.pdf")


    if make_movies:
        Q_movie = make_movie_imshow(Q_maps)
        sd_movie = make_movie_imshow(sd_maps)

        print("Writing movies (this can take some time and we cannot get progress)")
        Q_movie.save('Q_movie.mp4')
        sd_movie.save('sd_move.mp4')


if __name__ == "__main__":
    # Run in script mode

    import sys
    import os
    import pickle
    from multiprocessing import Pool
    from functools import partial

    # Physics Setup
    bbox_x = [-100, 100]
    bbox_y = bbox_x
    res_elem = 0.5
    res = survis.helper.get_res(res_elem, bbox_x, bbox_y)

    # Computing Setup
    n_cpus = os.cpu_count() if os.cpu_count() <= 8 else 8

    n_snaps = get_snaps()
    filenames = ["snapshot_{:03d}.hdf5".format(x) for x in range(n_snaps)]

    if "--test" in sys.argv:
        filenames = ['test_data.hdf5']

    if "--read" in sys.argv:
        print("Reading data")
        with open('processed_variables.pkl', 'rb') as pck:
            result = pickle.load(pck)

    else:
        print("Beginning data analysis \n")

        # We need more arguments than Pool.map allows
        t = tqdm(total=int(n_snaps/n_cpus) + 1, desc="Data Processing")

        def update():  # This is required because of pickling
            t.update()

        mapped_process = partial(processing_run,
                                 res=res, bbox_x=bbox_x, bbox_y=bbox_y,
                                 callback=update)

        with Pool(processes=n_cpus) as processing_pool:
            result = processing_pool.map(mapped_process, filenames)


    if "--save" in sys.argv:
        print("Saving data to processed_variables.pkl")
        with open('processed_variables.pkl', 'wb') as pck:
            pickle.dump(result, pck)


    if not ("--noplot" in sys.argv):
        print("Beginning data plotting")
        make_plots(result)
