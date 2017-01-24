# Preprocesses the gadget data into a useful format

import h5py
import numpy as np

def test():
    import matplotlib.pyplot as plt
    fname = "test_data.hdf5"
    res = 50
    DG = DataGridder(fname, res, res, -100, 100, -100, 100)

    print("Star Particle Attrs: {}".format([x for x in DG.star.keys()]))
    print("Gass Mass: {}; Star Mass: {}".format(DG.gas_mass, DG.star_mass))

    # Now we vis!
    gas_sd = np.reshape(DG.gas_mass_arr, (res, res))
    plt.imshow(gas_sd)
    plt.show()

    return


class DataGridder(object):
    def __init__(self, fname, binsx, binsy, xmin, xmax, ymin, ymax):
        # note that binsx and binsy should be similar to the smoothing
        # lengh used in the simulation.
        self.fname = fname

        self.binsx = binsx
        self.binsy = binsy
        self.xmin = xmin
        self.ymin = ymin
        self.xmax = xmax
        self.ymax = ymax

        self.header, self.gas, self.star = self.read_data()
        self.extract_header()
        self.gas_id_arr, self.gas_mass_arr = self.bin_data(self.gas, self.gas_mass)

        return


    def read_data(self):
        # Breaks down the data into particle types
        f = h5py.File(self.fname, 'r')

        # Header, gas, stars
        return f['Header'], f['PartType0'], f['PartType4']


    def extract_header(self):
        # Extracts the useful information in the Header
        self.time = self.header.attrs['Time']
        self.box_size = self.header.attrs['BoxSize']
        self.gas_mass = self.header.attrs['MassTable'][0]
        self.star_mass = self.header.attrs['MassTable'][4]

        return self.time, self.box_size, self.gas_mass, self.star_mass


    def bin_data(self, raw_data, part_mass):
        # raw_data is e.g. GADGET['PartType0'].
        # grids are left as flat lists for efficiency
        id_grid = [[] for x in range((self.binsx*self.binsy))]
        mass_arr = np.zeros(self.binsx*self.binsy)

        n_particles = len(raw_data['Coordinates'])

        coords = raw_data['Coordinates']
        ids = raw_data['ParticleIDs']

        binsize_x = (self.xmax - self.xmin)/(self.binsx)
        binsize_y = (self.ymax - self.ymin)/(self.binsy)


        for particle in range(n_particles):
            binx = int(np.floor((coords[particle][0] - self.xmin)/binsize_x))
            biny = int(np.floor((coords[particle][1] - self.ymin)/binsize_y))

            # todo: this should be done in a preprocessing loop
            if ((binx > 0) and (biny > 0) and (binx < self.binsx) and (biny < self.binsy)):
                this_bin = binx + self.binsx*biny

                # Now we do the processing for this particle
                id_grid[this_bin].append(ids[particle])
                mass_arr[this_bin] += part_mass
            else:
                pass  # The particle does not lie within the range

        return id_grid, mass_arr
