import survis.preprocess as pre
import survis.toomre as toom
import survis.helper as hp
import survis.fiducial as fid
import matplotlib.pyplot as plt
import numpy as np

print("Running preprocess.py tests")
# Plots a few things with the example data

fname = "test_data.hdf5"
res_elem = 5
res = hp.get_res(res_elem, [-100, 100], [-100, 100])
DG = pre.DataGridder(fname, res[0], res[1], -100, 100, -100, 100)

print("Gass Mass: {}; Star Mass: {}".format(DG.gas_mass, DG.star_mass))

def plot_circ(R):
    rs = np.ones(100) * R
    thetas = np.arange(0, 1, 0.01) * np.pi * 2
    x = rs * np.cos(thetas)
    y = rs * np.sin(thetas)
    return x, y

# Now we vis!
f = plt.figure(1)
gas_sd = DG.gas_data['masses']
plt.imshow(gas_sd)
plt.title('Gas mass')
xs, ys = plot_circ(13 * 40/100.)
plt.plot(xs + 20 , ys + 20)
xs, ys = plot_circ(17 * 40/100.)
plt.plot(xs + 20 , ys + 20)
f.show()

g = plt.figure(2)
gas_v = DG.gas_data['velocities']
plt.imshow(gas_v)
plt.title('Gas velocity')
g.show()

h = plt.figure(3)
gas_d = DG.gas_data['densities']
plt.imshow(gas_d)
plt.title('Gas density')
h.show()

print("Running fiducial.py tests")

sgas, sstar = fid.surface_density(DG, 15, 2)
print("Gas surface density: {}; Star: {}".format(sgas, sstar))
print("Gas Toomre Q: {}".format(fid.toomre_Q_gas(DG, 15, 2)))

print("Running toomre.py tests")

import matplotlib.cm as cm
i = plt.figure(4)

toomQ = hp.get_toomre_Q(DG, toom.sound_speed, res_elem)
cmap = cm.viridis
cmap.set_bad('white', -1.)

plt.imshow(toomQ, cmap)
plt.title('Gas Q')
plt.colorbar()
i.show()

input()  # keep figures alive
