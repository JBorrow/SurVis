# Preprocesses the gadget data into a useful format

import h5py
import numpy as np

def test():
    fname = "test_data.hdf5"
    DG = DataGridder(fname, 1, 2, 3, 4, 5, 6)
    head, gas, star = DG.read_data()

    print([x for x in star.keys()])
    print(gas['Density'][0])
    return


class DataGridder(object):
    def __init__(self, fname, binsx, binsy, xmin, xmax, ymin, ymax):
        self.fname = fname

        self.binsx = binsx
        self.binsy = binsy
        self.xmin = xmin
        self.ymin = ymin
        self.xmax = xmax
        self.ymax = ymax

        return


    def read_data(self):
        # Breaks down the data into particle types
        f = h5py.File(self.fname, 'r')

        # Header, gas, stars
        return f['Header'], f['PartType0'], f['PartType4']


    def bin_data(self, raw_data):
        # raw_data is e.g. GADGET['PartType0'].
        # grids are left as flat lists for efficiency
        self.id_grid = [[]*(self.binsx*self.binsy)]
        self.mass_arr = np.zeros(self.binsx*self.binsy)

        n_particles = len(raw_data['Coordinates'])

        coords = raw_data['Coordinates']
        dens = raw_data['Density']


        for particle in range(n_particles):
            binx = np.floor(coords[particle])


        return
