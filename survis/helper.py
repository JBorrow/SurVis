""" Contains helper functions for common plotting operations,
    such as finding the appropriate number of bins given a
    smoothing length. """

import survis.toomre as toom
import numpy as np


def get_res(smoothing, bbox_x, bbox_y):
    """ Smoothing should be a number, bbox_x = [minx, maxx] """

    def res(bbox, sm):
        return int((bbox_x[1] - bbox_x[0])/sm)

    return [res(bbox_x, smoothing), res(bbox_y, smoothing)]


def get_toomre_Q(DG, sound):
    """ DG is the datagrid, and sound is a callable function
        that returns the sound speed (see toomre.py).
        It automatically masks regions where there are no particles,
        which can be used through cmap.set_bad()."""

    gas_sd = DG.gas_data['masses']
    star_sd = DG.star_data['masses']

    gas_v = DG.gas_data['velocities']
    gas_d = DG.gas_data['densities']

    gas_q = toom.Q_gas(sound, gas_v, gas_d, gas_sd + star_sd)

    # now the sections with no particles must be masked

    return np.ma.array(gas_q, mask=(gas_q == 0.))
    


