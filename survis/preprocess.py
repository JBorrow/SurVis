"""
Preprocesses the gadget data into a useful format
Contains the object DataGridder which does the heavy lifting with
the function bin_data which iterates through the supplied file.
This data can then be visualised (see test() for an example).
"""

import h5py
import numpy as np


class DataGridder(object):
    def __init__(self, fname, binsx, binsy, xmin, xmax, ymin, ymax, autobin=True):
        """ note that binsx and binsy should be similar to the smoothing
            lengh used in the simulation.

            If you do not want the data to be autmatically binned on the
            initialization of DataGridder, set autobin to false."""

        self.fname = fname

        self.binsx = binsx
        self.binsy = binsy
        self.xmin = xmin
        self.ymin = ymin
        self.xmax = xmax
        self.ymax = ymax

        self.header, self.gas, self.star = self.read_data()
        self.extract_header()

        if autobin:
            self.gas_data = self.bin_data(self.gas, self.gas_mass, True)
            self.star_data = self.bin_data(self.star, self.star_mass, False)

        return


    def read_data(self):
        """ Breaks down the data into particle types """
        f = h5py.File(self.fname, 'r')

        # Header, gas, stars
        return f['Header'], f['PartType0'], f['PartType4']


    def extract_header(self):
        """ Extracts the useful information in the Header """
        self.time = self.header.attrs['Time']
        self.box_size = self.header.attrs['BoxSize']
        self.gas_mass = self.header.attrs['MassTable'][0]
        self.star_mass = self.header.attrs['MassTable'][4]

        return self.time, self.box_size, self.gas_mass, self.star_mass


    def bin_data(self, data, part_mass, hydro=True, ids=False):
        """ raw_data is e.g. GADGET['PartType0'].
            grids are left as flat lists for efficiency
            vel_grid returns v/r for each **particle** in a similar way to
            id_grid, use mean_grid to bin fully
            if(hyrdo) we also bin and return density (pressure)
            if (ids) we give a list of ids per bin (warning:SLOW) """

        vel_arr = np.zeros((self.binsx, self.binsy))
        n_arr = np.zeros((self.binsx, self.binsy))

        if (ids):
            id_grid = [[] for x in range((self.binsx*self.binsy))]

        if (hydro):
            d_arr = np.zeros((self.binsx, self.binsy))

        n_particles = len(data['Coordinates'])

        binsize_x = (self.xmax - self.xmin)/(self.binsx)
        binsize_y = (self.ymax - self.ymin)/(self.binsy)

        binsx = ((data['Coordinates'][:, 0] - self.xmin)/binsize_x).astype(int)
        binsy = ((data['Coordinates'][:, 1] - self.ymin)/binsize_y).astype(int)

        radii2 = np.sum(np.square(data['Coordinates'][()]), 1)
        vels = np.sqrt(np.sum(np.square(data['Velocities'][()]), 1)/radii2)

        if (ids):
            print("WARNING: The ID feature is not implemented")

        if(hydro):
            for v, bx, by, rho in zip(vels, binsx, binsy, data['Density']):
                try:
                    vel_arr[bx, by] += v
                    d_arr[bx, by] += rho
                    n_arr[bx, by] += 1
                except IndexError:
                    # Particle out of range, ignore
                    pass
        else:
            for v, bx, by in zip(vels, binsx, binsy):
                try:
                    vel_arr[bx, by] += v
                    n_arr[bx, by] += 1
                except IndexError:
                    # Particle out of range, ignore
                    pass

        # Tidy up

        m_arr = n_arr * part_mass

        # To prevent divide by 0 errors, we will have 0 velocity anyway
        n_arr[n_arr == 0] = 1
        vel_arr = np.sqrt(2) * vel_arr/n_arr

        if (hydro):
            d_arr = d_arr/n_arr

        ret = {'masses' : m_arr,
               'velocities' : vel_arr,}

        if (hydro):
            ret['densities'] = d_arr

        if (ids):
            ret['ids'] = id_grid

        return ret
