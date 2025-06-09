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

# Load the subhalo data from the TNG-300-2 simulation
gal_data = il.groupcat.loadSubhalos(basePatch, 99)

Flag = gal_data['SubhaloFlag']  #flag for valid subhalos
Mgal = gal_data['SubhaloMassType'][:, 1] * 10**10 / h  #Msun
x = gal_data['SubhaloPos'][:, 0] * 1e-3 / h  #Mpc
y = gal_data['SubhaloPos'][:, 1] * 1e-3 / h  #Mpc
z = gal_data['SubhaloPos'][:, 2] * 1e-3 / h  #Mpc
vx = gal_data['SubhaloVel'][:, 0]
vy = gal_data['SubhaloVel'][:, 1]
vz = gal_data['SubhaloVel'][:, 2]

# Filter out invalid subhalos
Mgal = Mgal[Flag == 1]  #only valid subhalos
x = x[Flag == 1]
y = y[Flag == 1]
z = z[Flag == 1]
vx = vx[Flag == 1]
vy = vy[Flag == 1]
vz = vz[Flag == 1]

#from 0 to L0, to -L0/2 to L0/2
x = x - L0/2
y = y - L0/2
z = z - L0/2

#######################################################################################################
#IN THIS CASE, vx,vy,vz are already in physical units, so no need to multiply by 1/a, as in haloes
#######################################################################################################

mcutdown = 1e7  #Msun
mcutup = 1e12  #Msun
mask = (Mgal > mcutdown) & (Mgal < mcutup)
npart = np.count_nonzero(mask)

#everything to int32, float32
x = x.astype(np.float32)[mask]
y = y.astype(np.float32)[mask]
z = z.astype(np.float32)[mask]
vx = vx.astype(np.float32)[mask]
vy = vy.astype(np.float32)[mask]
vz = vz.astype(np.float32)[mask]
Mgal = Mgal.astype(np.float32)[mask]
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
print('Mvir:', f'{np.min(Mgal):.2e}', f'{np.max(Mgal):.2e}')

#INPUT FILE FOR AVISM
#Write binary file with positions, velocities and mass
with FF('bin_file_part00001_TNG_galaxies', 'w') as f:
    f.write_record(npart, zeta) #header
    f.write_record(x)
    f.write_record(y)
    f.write_record(z)
    f.write_record(vx)
    f.write_record(vy)
    f.write_record(vz)
    f.write_record(Mgal)
