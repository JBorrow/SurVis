""" Contains objects used for analysis, primarily with common.py
    within this module. """

import numpy as np
import survis

class CommonDataObject(object):
    """ This object does a processing run and extracts a bunch of information
        from a given snapshot file, which is determined by filename. """

    def __init__(self, filename, res, bbox_x, bbox_y, elem_size):
        self.filename = filename

        self.res = res
        self.bbox_x = bbox_x
        self.bbox_y = bbox_y
        self.elem_size = elem_size

        # These may be modified later on before running analysis
        self.sound_speed = survis.toomre.sound_speed_sne
        self.solar_radius = 8
        self.smoothing = 0.4

        return


    def run_analysis(self):
        data_grid = survis.preprocess.DataGridder(self.filename,
                                                  self.res[0],
                                                  self.res[1],
                                                  self.bbox_x[0],
                                                  self.bbox_x[1],
                                                  self.bbox_y[0],
                                                  self.bbox_y[1])

        self.Q_map = survis.helper.get_toomre_Q(data_grid,
                                                self.sound_speed,
                                                self.elem_size)

        # Normally the masses of each element are given, we must divide by size
        # as well as a conversion factor to give Msun / pc^2
        self.sd_map = data_grid.gas_data['masses']/((1e6) * self.elem_size**2)


        # Now the values at a given radius
        self.sd_r = survis.fiducial.surface_density(data_grid,
                                                    self.solar_radius,
                                                    self.smoothing)
        self.Q_r = survis.fiducial.toomre_Q_gas(data_grid,
                                                self.solar_radius,
                                                self.smoothing,
                                                self.sound_speed)

        # Now the values for all radii
        self.Q_variation_with_r = survis.helper.toomre_Q_r(data_grid,
                                                           self.sound_speed,
                                                           self.smoothing,
                                                           self.bbox_x[1])

        self.sd_variation_with_r = survis.helper.sd_r(data_grid,
                                                      self.smoothing,
                                                      self.bbox_x[1])

        self.n_part_r, self.bins = survis.helper.n_particles_bins(data_grid)

        self.vert_opt, self.vert_err = survis.profiles.vertical_profile(data_grid)

        return


class CommonDataExtractor(object):
    """ This object is used to extract the data back to lists per snapshot.
        We begin with [CommonDataObject, CommonDataObject, ...] but really
        we want the actual data items [snap0, snap1, snap2] x N. This does
        that. """

    def __init__(self, cdo_list):
        self.cdo_list = cdo_list

        self.Q_map = []
        self.sd_map = []
        self.sd_r = []
        self.Q_r = []
        self.Q_variation_with_r = []
        self.sd_variation_with_r = []
        self.n_part_r = []
        self.bins = []
        self.vert_opt = []
        self.vert_err = []

        self._reshape()

        self.sd_r = self.clean(self.sd_r)
        self.Q_r = self.clean(self.Q_r)
        self.n_part_r = self.clean(self.n_part_r)

        return

    
    def clean(self, my_list):
        """ Some lists turn up in a bad shape. This cleans them """
        return np.array([np.array(x) for x in my_list])

    
    def _reshape(self):
        for item in tqdm(self.cdo_list, desc="Reshaping data"):
            self.Q_map.append(item.Q_map)
            self.sd_map.append(item.sd_map)
            self.sd_r.append(item.sd_r)
            self.Q_r.append(item.Q_r)
            self.Q_variation_with_r.append(item.Q_variation_with_r)
            self.sd_variation_with_r.append(item.sd_variation_with_r)
            self.n_part_r.append(item.n_part_r)
            self.bins.append(item.bins)
            self.vert_opt.append(item.vert_opt)
            self.vert_err.append(item.vert_err)

        return


