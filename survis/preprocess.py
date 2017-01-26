# Preprocesses the gadget data into a useful format
# Contains the object DataGridder which does the heavy lifting with
# the function bin_data which iterates through the supplied file.
# This data can then be visualised (see test() for an example).

import h5py
import numpy as np


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

        self.gas_data = self.bin_data(self.gas, self.gas_mass, True)
        self.star_data = self.bin_data(self.star, self.star_mass, False)

        return


    def read_data(self):
        # Breaks down the data into particle types
        f = h5py.File(self.fname, 'r')

        # Header, gas, stars
        return f['Header'], self.clean_data(f['PartType0'], True), self.clean_data(f['PartType4'])


    def clean_data(self, data, hydro=False):
        # retrns the gadget data in an iteratable *BY PARTICLE* format
        n_particles = len(data['Coordinates'])

        if (hydro):
            ret = np.array([data['Coordinates'],
                            data['Velocities'],
                            data['ParticleIDs'],
                            data['Density']])

        else:
            ret = np.array([data['Coordinates'],
                            data['Velocities'],
                            data['ParticleIDs']])

        return ret.T


    def extract_header(self):
        # Extracts the useful information in the Header
        self.time = self.header.attrs['Time']
        self.box_size = self.header.attrs['BoxSize']
        self.gas_mass = self.header.attrs['MassTable'][0]
        self.star_mass = self.header.attrs['MassTable'][4]

        return self.time, self.box_size, self.gas_mass, self.star_mass


    def radii(self, coords):
        # Takes the coordinates and returns the radii
        return np.sqrt(coords[0]**2 + coords[1]**2 + coords[2]**2)


    def rms(self, item):
        return np.sqrt(np.mean(np.square(item)))


    def mean_grid(self, grid):
        # Takes a grid (well, a flattened one, such as id_grid) and means
        # it across the list
        return [(sum(x)/len(x)) if len(x) > 0 else 0. for x in grid]


    def bin_data(self, raw_data, part_mass, hydro=True, ids=False):
        # raw_data is e.g. GADGET['PartType0'].
        # grids are left as flat lists for efficiency
        # vel_grid returns v/r for each **particle** in a similar way to
        # id_grid, use mean_grid to bin fully
        # if(hyrdo) we also bin and return density (pressure)
        # if (ids) we give a list of ids per bin (warning:SLOW)

        vel_arr = np.zeros(self.binsx*self.binsy)
        n_arr = np.zeros(self.binsx*self.binsy)

        if (ids):
            id_grid = [[] for x in range((self.binsx*self.binsy))]

        if (hydro):
            d_arr = np.zeros(self.binsx*self.binsy)

        n_particles = len(raw_data)

        binsize_x = (self.xmax - self.xmin)/(self.binsx)
        binsize_y = (self.ymax - self.ymin)/(self.binsy)

        for particle in raw_data:
            binx = int(np.floor((particle[0][0] - self.xmin)/binsize_x))
            biny = int(np.floor((particle[0][1] - self.ymin)/binsize_y))

            # todo: this should really be done in a preprocessing loop
            if ((binx > 0) and (biny > 0) and (binx < self.binsx) and (biny < self.binsy)):
                this_bin = binx + self.binsx*biny

                n_arr[this_bin] += 1
                vel_arr[this_bin] += (self.rms(particle[1])/self.radii(particle[1]))/(3.086e16)
                
                if (ids):
                    id_grid[this_bin].append(particle[2])

                if (hydro):
                    d_arr[this_bin] += particle[3]

            else:
                pass  # The particle does not lie within the range

        # Tidy up
        
        m_arr = n_arr * part_mass

        # To prevent divide by 0 errors, we will have 0 velocity anyway
        n_arr[n_arr == 0] = 1
        vel_arr = vel_arr/n_arr
        
        if (hydro):
            d_arr = d_arr/n_arr

        
        ret = {'masses' : m_arr,
               'velocities' : vel_arr,}

        if (hydro):
            ret['densities'] = d_arr

        if (ids):
            ret['ids'] = id_grid

        return ret


