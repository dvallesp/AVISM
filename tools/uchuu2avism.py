"""
Example script to handle particle simulation data and write a binary file (particles) to be read by AVISM.
In this particular case, dark matter haloes were utilised.
The extrapolation to raw particle simulation data is straightforward.
"""

import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
from scipy.io import FortranFile as FF

h = 0.6774 #uchuu
L0 = 400./h
c = 299792.458
zeta = 0.00 #redshift

#path to uchuu_files
path = 'uchuu_files'

file = 'MiniUchuu_halolist_z0p00.h5'
file = os.path.join(path, file)

# Open the H5 file in read mode
with h5py.File(file, 'r') as file:
    Mvir = file['Mvir'][:]
    x = file['x'][:]
    y = file['y'][:]
    z = file['z'][:]
    vx = file['vx'][:]
    vy = file['vy'][:]
    vz = file['vz'][:]

#from Mpc/h to Mpc
x = x/h
y = y/h
z = z/h

#from 0 to L0, to -L0/2 to L0/2
x = x - L0/2
y = y - L0/2
z = z - L0/2

#Cut just for Mvir>10^10
mask = Mvir > 1e10
x = x[mask]
y = y[mask]
z = z[mask]
vx = vx[mask]
vy = vy[mask]
vz = vz[mask]
Mvir = Mvir[mask]

npart = len(x)


#everything to int32, float32
x = x.astype(np.float32)
y = y.astype(np.float32)
z = z.astype(np.float32)
vx = vx.astype(np.float32)
vy = vy.astype(np.float32)
vz = vz.astype(np.float32)
Mvir = Mvir.astype(np.float32)
npart = np.int32(npart)

#CHECK!!
print('Checking data...')
print('L0', L0)
print('npart:', npart)
print('x:', np.min(x), np.max(x))
print('y:', np.min(y), np.max(y))
print('z:', np.min(z), np.max(z))
print('vx:', np.min(vx), np.max(vx))
print('vy:', np.min(vy), np.max(vy))
print('vz:', np.min(vz), np.max(vz))
print('Mvir:', f'{np.min(Mvir):.2e}', f'{np.max(Mvir):.2e}')

#INPUT FILE FOR AVISM
#Write binary file with positions, velocities and mass
with FF('bin_file_part00001', 'w') as f:
    f.write_record(npart, zeta) #header
    f.write_record(x)
    f.write_record(y)
    f.write_record(z)
    f.write_record(vx)
    f.write_record(vy)
    f.write_record(vz)
    f.write_record(Mvir)

#visual inspection
# fig, ax = plt.subplots(1, 1, figsize=(5, 4), dpi = 200)
# mask1 = z < 1
# mask2 = z > -1
# mask = mask1 & mask2
# ax.scatter(x[mask], y[mask], s=0.1)
# plt.savefig('part.png')