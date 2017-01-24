import survis.preprocess as pre
import matplotlib.pyplot as plt
import numpy as np

print("Running preprocess.py tests")
# Plots a few things with the example data

fname = "test_data.hdf5"
res = 50
DG = pre.DataGridder(fname, res, res, -100, 100, -100, 100)

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

h = plt.figure(3)
gas_d = np.reshape(DG.mean_gas_d_arr, (res, res))
plt.imshow(gas_d)
plt.title('Gas density')
h.show()
input()  # keep figures alive


