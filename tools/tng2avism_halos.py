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
import illustris_python as il

#FOR A DESCRIPTION OF ILLUSTRIS-TNG VARIABLES:
# https://www.illustris-project.org/data/docs/specifications/

h = 0.6774 
L0 = 302.63
c = 299792.458
zeta = 0.00 #redshift

basePatch = '/scratch/monllor/TNG-300-2'
fields = ['']

# Load the halo data from the TNG-300-2 simulation
halo_data = il.groupcat.loadHalos(basePatch, 99)

Mvir = halo_data['Group_M_Crit200'] * 10**10 / h  #Msun
x = halo_data['GroupPos'][:, 0] * 1e-3 / h  
y = halo_data['GroupPos'][:, 1] * 1e-3 / h  
z = halo_data['GroupPos'][:, 2] * 1e-3 / h  
vx = halo_data['GroupVel'][:, 0]
vy = halo_data['GroupVel'][:, 1] 
vz = halo_data['GroupVel'][:, 2]


#from 0 to L0, to -L0/2 to L0/2
x = x - L0/2
y = y - L0/2
z = z - L0/2

#################################################################################
#TAKE INTO ACCOUNT vx,vy,vz should be multiplied by 1/a, but since this is z = 0, a = 1
#################################################################################

mcut = 1e9  #Msun
npart = np.count_nonzero(Mvir > mcut)

#everything to int32, float32
x = x.astype(np.float32)[Mvir > mcut]
y = y.astype(np.float32)[Mvir > mcut]
z = z.astype(np.float32)[Mvir > mcut]
vx = vx.astype(np.float32)[Mvir > mcut]
vy = vy.astype(np.float32)[Mvir > mcut]
vz = vz.astype(np.float32)[Mvir > mcut]
Mvir = Mvir.astype(np.float32)[Mvir > mcut]
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
with FF('bin_file_part00001_TNG_halos', 'w') as f:
    f.write_record(npart, zeta) #header
    f.write_record(x)
    f.write_record(y)
    f.write_record(z)
    f.write_record(vx)
    f.write_record(vy)
    f.write_record(vz)
    f.write_record(Mvir)