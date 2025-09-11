"""
Script to handle Manticore-Local data and write a binary file to be read by AVISM.
Manticore-Local: https://arxiv.org/abs/2505.10682
"""

import numpy as np
import os
import h5py
from scipy.io import FortranFile as FF

h = 0.681
L0 = 1000. #cMpc
nx = 256
dx = 3.9 #cMpc
zeta = 0.0 #redshift
R_mask = 200/h
path = '/scratch/monllor/Manticore-Local/N256'

#arrays of realisations
nsim = 80
dens_array = np.zeros((nsim, nx, nx, nx))
vx_array = np.zeros((nsim, nx, nx, nx))
vy_array = np.zeros((nsim, nx, nx, nx))
vz_array = np.zeros((nsim, nx, nx, nx))

for i in range(nsim):
    #Open hdf5 files
    #Keys: 'density', 'num_in_cell', 'p0', 'p1', 'p2'
    filename = f'mcmc_{i}.hdf5'
    file = os.path.join(path, filename)
    with h5py.File(file, 'r') as f:
        density = f['density'][:]
        p0 = f['p0'][:]
        p1 = f['p1'][:]
        p2 = f['p2'][:]
        vx = p0 / density
        vy = p1 / density
        vz = p2 / density
        mean_dens = np.mean(density)
        dens = density/mean_dens

        dens_array[i] = dens
        vx_array[i] = vx
        vy_array[i] = vy
        vz_array[i] = vz    

#Mean over realisations
dens_mean = np.mean(dens_array, axis=0)
vx_mean = np.mean(vx_array, axis=0)
vy_mean = np.mean(vy_array, axis=0)
vz_mean = np.mean(vz_array, axis=0)
delta_mean = dens_mean - 1.


#Outside R = R_mask, put zeros
mask = np.ones((nx, nx, nx), dtype=bool)
for i in range(nx):
    for j in range(nx):
        for k in range(nx):
            x = (i - nx/2)*dx
            y = (j - nx/2)*dx
            z = (k - nx/2)*dx
            r = np.sqrt(x**2 + y**2 + z**2)
            if r > R_mask:
                mask[i,j,k] = False

delta_mean = np.where(mask, delta_mean, 0.)

#CHECK!
#vel in km/s
print('delta', np.max(delta_mean), np.min(delta_mean))
print('vx', np.max(vx_mean), np.min(vx_mean))
print('vy', np.max(vy_mean), np.min(vy_mean))
print('vz', np.max(vz_mean), np.min(vz_mean))

with FF('bin_file_grid00001', 'w') as f:
    f.write_record(np.float32(zeta))
    f.write_record(np.float32(delta_mean.T))
    f.write_record(np.float32(vx_mean.T))
    f.write_record(np.float32(vy_mean.T))
    f.write_record(np.float32(vz_mean.T))