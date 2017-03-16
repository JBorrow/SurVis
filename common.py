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

def processing_run(filename, res, bbox_x, bbox_y, elem_size, callback=None):
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
                                       elem_size)

    # Normally the masses of each element are given, we must divide by size
    # as well as a conversion factor to give Msun / pc^2
    sd_map = data_grid.gas_data['masses']/((1e6) * elem_size**2)


    # Now the values at a given radius
    sd_r = survis.fiducial.surface_density(data_grid, solar_radius, smoothing)
    Q_r = survis.fiducial.toomre_Q_gas(data_grid,
                                       solar_radius,
                                       smoothing,
                                       survis.toomre.sound_speed_sne)

                                       # Now the values for all radii
    Q_variation_with_r = survis.helper.toomre_Q_r(data_grid,
                                                  survis.toomre.sound_speed_sne,
                                                  smoothing,
                                                  bbox_x[1])

    n_part_r, bins = survis.helper.n_particles_bins(data_grid)

    if not (callback is None):
        callback()

    return Q_map, sd_map, Q_r, sd_r, Q_variation_with_r, n_part_r, bins, filename[9:12]


def get_snaps(directory = "."):
    n_snaps = 0

    while "snapshot_{:03d}.hdf5".format(n_snaps) in os.listdir(directory):
        n_snaps += 1
        if n_snaps > 10000:
            break

    return n_snaps


def make_movie_imshow(data, bad_color='black', log=False, vmin=0, vmax=3):
    images = []

    fig = plt.figure()
    colormap = cm.get_cmap('viridis')
    c_scale = col.Normalize(vmin=vmin , vmax=vmax)
    colormap.set_bad(bad_color, 1.0)


    for item in tqdm(data, desc="Data plotting"):
        images.append([plt.imshow(item, cmap=colormap, vmin=vmin, vmax=vmax)])

    fig.colorbar(cmap=colormap, norm=c_scale, mappable=images[-1][0])

    return animation.ArtistAnimation(fig, images, interval=50, repeat_delay=3000,
                                     blit=True)


def make_linear_plot_movie(data, ylabel, ymin=0, ymax=0):
    xs = np.arange(len(data[0]))*smoothing
    n_images = len(data)
    images = []

    fig, ax = plt.subplots()
    ax.plot([0, 1000], [1, 1], 'k--')

    for item in tqdm(data, desc="Linear plot"):
        images.append(ax.plot(xs, item, 'b-'))

    plt.xlim([0, 100*smoothing])
    plt.ylim([0, 1.5])

    return animation.ArtistAnimation(fig, images, interval=50, repeat_delay=3000,
                                     blit=True)


def make_linear_plot(data, ylabel, ymin=0, ymax=5.):
    n_snaps = len(data)
    fig, ax = plt.subplots()

    ax.plot(np.arange(n_snaps), np.array(data))

    ax.set_xlim(0, n_snaps)
    #ax.set_ylim(ymin, ymax)

    ax.set_xlabel("Snapshot number")
    ax.set_ylabel(ylabel)

    return fig, ax


def n_part_r_plot(n_r, bin_edges):
    n_snaps = len(n_r[0, :])
    fig, ax = plt.subplots()

    for index, dataset in enumerate(n_r):
        label = "{} $\leq r <$ {}".format(bin_edges[i], bin_edges[i+1]) 
        ax.plot(range(n_snaps), dataset, label=label)

    ax.set_xlabel("Snapshot number")
    ax.set_xlims(0, n_snaps)
    ax.set_ylabel("Number of particles within bounds")

    ax.legend()

    return fig, ax



def make_plots(result, make_movies=True, show_plots=False):
    result = np.array(result)
    Q_maps = result.T[0]
    sd_maps = result.T[1]
    Q_r = result.T[2]
    sd_r = np.array([np.array(x) for x in result.T[3]]).T
    Q_variation_with_r = result.T[4]
    n_part_r = result.T[5]
    bin_edges = result.T[6]
    # We have to manually convert the lists here from when they get pickled
    sd_r_gas = sd_r[0, :]
    sd_r_star = sd_r[1, :]

    Q_fig, Q_ax = make_linear_plot(Q_r, "Toomre $Q$", 0.5, 3.0)
    sd_fig, sd_ax = make_linear_plot(sd_r_gas, "Surface Density ($M_\odot$ pc$^{-2}$)", 0, 1e5)
    n_fig, n_ax = n_part_r_plot(n_part_r, bin_edges)

    if show_plots:
        Q_fig.show()
        sd_fig.show()
        n_fig.show()
        print("Showing plots, press any key to continue")
        input()
    else:
        print("Writing plots")
        Q_fig.savefig("Q_fig.pdf")
        sd_fig.savefig("sd_gas.pdf")
        n_fig.savefig("n_gas.pdf")


    if make_movies:
        Q_movie = make_movie_imshow(Q_maps, vmin=0, vmax=2)
        sd_movie = make_movie_imshow(sd_maps, vmin=0, vmax=50)
        Q_of_r_mov = make_linear_plot_movie(Q_variation_with_r, "Q")

        print("Writing movies (this can take some time and we cannot get progress)")
        Q_movie.save('Q_movie.mp4')
        sd_movie.save('sd_movie.mp4')
        Q_of_r_mov.save('Q_of_r_mov.mp4')


if __name__ == "__main__":
    # Run in script mode

    import sys
    import os
    import pickle
    from multiprocessing import Pool
    from functools import partial

    # Physics Setup
    bbox_x = [-30, 30]
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
                                 elem_size=res_elem,
                                 callback=update)

        with Pool(processes=n_cpus) as processing_pool:
            result = processing_pool.map(mapped_process, filenames)


    if "--save" in sys.argv:
        print("Saving data to processed_variables.pkl")
        with open('processed_variables.pkl', 'wb') as pck:
            pickle.dump(result, pck)

    if not ("--noplot" in sys.argv):
        print("Beginning data plotting")
        show_plots = "--showplots" in sys.argv
        make_plots(result, show_plots=show_plots)
