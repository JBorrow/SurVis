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
    f = plt.figure(1)
    gas_sd = np.reshape(DG.gas_mass_arr, (res, res))
    plt.imshow(gas_sd)
    plt.title('Gas mass')
    f.show()

    g = plt.figure(2)
    gas_v = np.reshape(DG.mean_gas_vel_arr, (res, res))
    plt.imshow(gas_v)
    plt.title('Gas velocity')
    g.show()

    input()  # keep figures alive

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
        self.gas_id_arr, self.gas_vel_arr, self.gas_mass_arr = self.bin_data(self.gas, self.gas_mass)
        self.mean_gas_vel_arr = self.mean_grid(self.gas_vel_arr)

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


    def radii(self, coords):
        # Takes the coordinates and returns the radii
        return np.sqrt(coords[:, 0]**2 + coords[:, 1]**2 + coords[:, 2]**2)


    def rms(self, item):
        return np.sqrt(np.mean(np.square(item)))


    def mean_grid(self, grid):
        # Takes a grid (well, a flattened one, such as id_grid) and means
        # it across the list
        return [(sum(x)/len(x)) if len(x) > 0 else 0. for x in grid]


    def bin_data(self, raw_data, part_mass):
        # raw_data is e.g. GADGET['PartType0'].
        # grids are left as flat lists for efficiency
        # vel_grid returns v/r for each **particle** in a similar way to
        # id_grid, use mean_grid to bin fully

        id_grid = [[] for x in range((self.binsx*self.binsy))]
        vel_grid = [[] for x in range((self.binsx*self.binsy))]
        mass_arr = np.zeros(self.binsx*self.binsy)

        n_particles = len(raw_data['Coordinates'])

        coords = raw_data['Coordinates']
        ids = raw_data['ParticleIDs']
        vels = raw_data['Velocities']
        radiis = self.radii(coords)

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
                vel_grid[this_bin].append(self.rms(vels[particle])/radiis[particle])
                mass_arr[this_bin] += part_mass
            else:
                pass  # The particle does not lie within the range

        return id_grid, vel_grid, mass_arr
