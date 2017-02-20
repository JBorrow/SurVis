import survis.preprocess as pre
import survis.toomre as toom
import survis.helper as hp
import matplotlib.pyplot as plt
import numpy as np

print("Running preprocess.py tests")
# Plots a few things with the example data

fname = "test_data.hdf5"
res = 50
DG = pre.DataGridder(fname, res, res, -100, 100, -100, 100)

print("Gass Mass: {}; Star Mass: {}".format(DG.gas_mass, DG.star_mass))

# Now we vis!
f = plt.figure(1)
gas_sd = DG.gas_data['masses']
plt.imshow(gas_sd)
plt.title('Gas mass')
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

print("Running toomre.py tests")

import matplotlib.cm as cm
i = plt.figure(4)

toomQ = hp.get_toomre_Q(DG, toom.sound_speed)
cmap = cm.viridis
cmap.set_bad('white', -1.)

plt.imshow(toomQ, cmap)
plt.title('Gas Q')
plt.colorbar()
i.show()

input()  # keep figures alive
