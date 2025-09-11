"""
Script to handle Carrick et al. (2015)  data and write a binary file to be read by AVISM.
Pub: Carrick, J., Turnbull, S. J., Lavaux, G., & Hudson, M. J. 2015, Monthly Notices of the Royal Astronomical Society, 450, 317
"""

import numpy as np
import os
from scipy.io import FortranFile as FF

path = '2M++-linear-input'

# files
densfile = 'twompp_density.npy'
velfile  = 'twompp_velocity.npy'

densfile = os.path.join(path, densfile)
velfile = os.path.join(path, velfile)

# parameters
L = 587.3715 #Mpc, 400 Mpc/h  (h = 0.681)
nx = 257
ny = 257
nz = 257
c = 299792.458 #km/s
dx = L/nx
zeta = 0.

delta = np.float32(np.load(densfile))
vel = np.float32(np.load(velfile))
velx = vel[0,:,:,:]
vely = vel[1,:,:,:]
velz = vel[2,:,:,:]

#CHECK!!
#vel in km/s
print('Checking data...')
print('delta: ', np.max(delta), np.min(delta))
print('velx: ', np.max(velx), np.min(velx))
print('vely: ', np.max(vely), np.min(vely))
print('velz: ', np.max(velz), np.min(velz))


#INPUT FILE FOR AVISM
#Write binary file with delta a velocities
with FF('bin_file_grid00001', 'w') as f:
    f.write_record(np.float32(zeta))
    f.write_record(delta.T)
    f.write_record(velx.T)
    f.write_record(vely.T)
    f.write_record(velz.T)
