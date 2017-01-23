# Preprocesses the gadget data into a useful format

import h5py
import numpy as np

def test():
    fname = "test_data.hdf5"
    head, gas, star = read_data(fname)

    print([x for x in star.keys()])
    print(gas['Density'][0])
    return


def read_data(fname):
    # Breaks down the data into particle types
    f = h5py.File(fname, 'r')

    # Header, gas, stars
    return f['Header'], f['PartType0'], f['PartType4']


def bin_data(raw_data, binsx, binsy, xmin, xmax, ymin, ymax):
    # raw_data is e.g. GADGET['PartType0'].
    # grids are left as flat lists for efficiency
    id_grid = [[]*(binsx*binsy)]
    mass_arr = np.zeros(binsx*binsy)

    n_particles = len(raw_data['Coordinates'])

    coords = raw_data['Coordinates']
    dens = raw_data['Density']


    for particle in range(n_particles):
        binx = np.floor(coords[particle])


    return
