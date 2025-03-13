"""
Script to handle galaxy survey data and write a binary file (particles) to be read by AVISM.
Thought for all-sky surveys.
Steps:
- Read the galaxy survey data
- Translate to cartesian coordinates
- Perform a cut in redshift space (radial cut)
- Estimate the mass of the tracers calculating mass bias by means of background density
- Write a binary file to be read by AVISM
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from numba import njit
from scipy.io import FortranFile as FF
from scipy.ndimage import gaussian_filter1d

#astropy
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM, z_at_value
import astropy.units as u

#input option
type_file = 0 #0 for fits

#path to the galaxy survey data
survey_path = 'survey_input'
name = '2mrs_1175_done'

#constants and cosmology
c = 299792. #km/s
h = 0.678
omega_m = 0.31
cosmo = FlatLambdaCDM(H0=100*h, Om0=omega_m)


#read data
#####################
if type_file == 0:
#####################
    fitsfile = name+'.fits'
    fitsfile = os.path.join(survey_path, fitsfile)
    hdu = fits.open(fitsfile, dtype=np.float32)
    #BinTableHDU
    hdu = hdu[1]

#####################

#print keys
print('Keys : ')
print(hdu.columns.names)
print()

data = hdu.data
ra = data['RA']
dec = data['DEC']
vrec = data['V']
zeta = vrec / c
#########################################


#####################################
#comoving distance for each galaxy
#####################################
dc = cosmo.comoving_distance(zeta)
dc = dc.value

###############################
# n(R) profile
###############################
nbins = 100
radii_edge = np.linspace(0., 1.1*np.max(dc), nbins+1)
radii = 0.5*(radii_edge[1:] + radii_edge[:-1])
dr = radii[1] - radii[0]
counter = np.histogram(dc, bins=radii_edge)[0]
n =  counter / (4.*np.pi*radii**2*dr)

#plot
fig, axs = plt.subplots(1, 2, figsize=(8, 6), dpi = 200)
ax = axs[0]
ax.plot(radii, n, marker='o', color='black')
ax.set_xlabel('Radius (Mpc)')
ax.set_ylabel('Number density (1/Mpc^3)')
ax.set_yscale('log')
ax.set_xscale('log')
ax= axs[1]
ax.plot(radii, counter, marker='o', color='black')
ax.set_xlabel('Radius (Mpc)')
ax.set_ylabel('Number of galaxies')
ax.set_xscale('log')
ax.set_yscale('log')

plt.tight_layout()
plt.savefig('rho_survey.png')
##############################

###############################
#perform a cut in redshift space
z_cut = 0.068
d_cut = cosmo.comoving_distance(z_cut).value

print('Cutting at z : ', z_cut, ' d (Mpc) : ', d_cut)
print('Number of galaxies before cut : ', len(ra))

ra = ra[zeta < z_cut]
dec = dec[zeta < z_cut]
dc = dc[zeta < z_cut]
ngal = len(ra)

print('Number of galaxies after cut : ', ngal)
###############################

###############################
#convert to cartesian coordinates
coords = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, distance=dc)
###############################

#extract x, y, z
x, y, z = coords.cartesian.x.value, coords.cartesian.y.value, coords.cartesian.z.value

print()
print()

#PLOTS
#3D distribution
fig = plt.figure(figsize=(8, 8), dpi = 200)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, s=0.1)
ax.grid(False)
plt.savefig('galaxy_survey.png')


#####################################
#Tracer mass estimation using n(R) profile
#####################################

z_mean = np.nanmean(zeta)

#matter background density at z_mean
rho_mean = omega_m*cosmo.critical_density(z_mean).to(u.solMass/u.Mpc**3).value
print('Background density (Msun/Mpc^3): ', f'{rho_mean:.2e}')
Mtracer = np.zeros(ngal, dtype=np.float32)

for itracer in range(ngal):
    binn = np.argmin(np.abs(radii - dc[itracer]))
    M = rho_mean/n[binn]
    Mtracer[itracer] = M

print('Tracer mass range (Msun) : ', f'{np.min(Mtracer):.2e}', f'{np.max(Mtracer):.2e}')
print()
print()

####################################
#Border ghost particles
# - Place ghost particles outside the d_cut sphere
# - Mass is much larger than the tracers to avoid voids outside
#####################################

print('Generating ghost particles...')
nghost = 100_000
Mghost = 10**20 * np.ones(nghost, dtype=np.float64)

print('Ghost particles mass (Msun) : ', f'{Mghost[0]:.2e}')
print('Number of ghost particles : ', nghost)
print()
print()

def ghost_particles(nghost, d_cut):
    """
    Generate ghost particles inside the 2*d_cut cube
    But outside the d_cut sphere
    """
    xghost = np.zeros(nghost, dtype=np.float64)
    yghost = np.zeros(nghost, dtype=np.float64)
    zghost = np.zeros(nghost, dtype=np.float64)
    for i in range(nghost):
        xghost[i] = np.random.uniform(-d_cut, d_cut)
        yghost[i] = np.random.uniform(-d_cut, d_cut)
        zghost[i] = np.random.uniform(-d_cut, d_cut)
        while np.sqrt(xghost[i]**2 + yghost[i]**2 + zghost[i]**2) <= d_cut:
            xghost[i] = np.random.uniform(-d_cut, d_cut)
            yghost[i] = np.random.uniform(-d_cut, d_cut)
            zghost[i] = np.random.uniform(-d_cut, d_cut)
    return xghost, yghost, zghost

xghost, yghost, zghost = ghost_particles(nghost, d_cut)

### Concatenate tracers and ghost particles
#Types are float32 and int64 for npart
zeta = np.float32(z_mean)
npart = ngal + nghost
x = np.concatenate((x, xghost)).astype(np.float32)
y = np.concatenate((y, yghost)).astype(np.float32)
z = np.concatenate((z, zghost)).astype(np.float32)
Mtracer = np.concatenate((Mtracer, Mghost)).astype(np.float32)

#############
# SAVE TO BINARY FORTRAN FILE
#############
#Write binary file with positions, velocities and mass
with FF('bin_file_part00001', 'w') as f:
    f.write_record(npart, zeta) #header
    f.write_record(x)
    f.write_record(y)
    f.write_record(z)
    f.write_record(Mtracer)





