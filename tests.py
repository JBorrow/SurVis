import survis.preprocess as pre
import survis.toomre as toom
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
gas_sd = np.reshape(DG.gas_data['masses'], (res, res))
plt.imshow(gas_sd)
plt.title('Gas mass')
f.show()

g = plt.figure(2)
gas_v = np.reshape(DG.gas_data['velocities'], (res, res))
plt.imshow(gas_v)
plt.title('Gas velocity')
g.show()

h = plt.figure(3)
gas_d = np.reshape(DG.gas_data['densities'], (res, res))
plt.imshow(gas_d)
plt.title('Gas density')
h.show()

print("Running toomre.py tests")

i = plt.figure(4)
star_sd = np.reshape(DG.star_data['masses'], (res, res))
gas_Q = toom.Q_gas(toom.sound_speed, gas_v, gas_d, gas_sd + star_sd)
plt.imshow(gas_Q)
plt.title('Gas Q')
plt.colorbar()
i.show()

input()  # keep figures alive
