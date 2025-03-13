"""
Script to handle CORAS data and write a binary file to be read by AVISM.
CORAS: https://github.com/rlilow/CORAS
Pub: Robert Lilow & Adi Nusser, MNRAS 507, 1557–1581 (2021)
"""

import numpy as np
import os
from scipy.io import FortranFile as FF
import sys

#path to the CORAS data
coras_path = 'coras_input'

#CORAS files
densfile = 'cartesian_grid_density_zCMB.dat'
velfile  = 'cartesian_grid_velocity_zCMB.dat'

densfile = os.path.join(coras_path, densfile)
velfile = os.path.join(coras_path, velfile)

# CORAS parameters
L = 590. #Mpc, 400 Mpc/h  (h = 0.678)
nx = 201
ny = 201
nz = 201
sigma8 = 0.81
c = 299792.458 #km/s
dx = L/nx
zeta = 0.

#
ncells = nx*ny*nz

#READING CORAS FILES
#density
aux = np.zeros((ncells), dtype=np.float32)
with open(densfile, 'rb') as f:
    f.readline() #header
    for i in range(ncells):
        aux[i] = float(f.readline()[:-1])

delta = np.zeros((nx,ny,nz), dtype=np.float32)
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            delta[i,j,k] = aux[(i * nx + j) * ny + k]

delta = delta*sigma8

#velocity
auxx = np.zeros((ncells), dtype=np.float32)
auxy = np.zeros((ncells), dtype=np.float32)
auxz = np.zeros((ncells), dtype=np.float32)
with open(velfile, 'rb') as f:
    f.readline() #header
    for i in range(ncells):
        line = f.readline()
        #line are a bytes list of 3 floats separated by \t:
        numbers = line.split(b'\t')
        auxx[i] = float(numbers[0])
        auxy[i] = float(numbers[1])
        auxz[i] = float(numbers[2])

velx = np.zeros((nx,ny,nz), dtype=np.float32)
vely = np.zeros((nx,ny,nz), dtype=np.float32)
velz = np.zeros((nx,ny,nz), dtype=np.float32)

for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            velx[i,j,k] = auxx[(i * nx + j) * ny + k]
            vely[i,j,k] = auxy[(i * nx + j) * ny + k]
            velz[i,j,k] = auxz[(i * nx + j) * ny + k]


#CHECK!!
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
