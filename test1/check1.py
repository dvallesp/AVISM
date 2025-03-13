"""
THIS CODE IS FOR ASSESING THE ACCURACY OF THE VOID FINDER ALGORITHM
IN FINDING VOIDS IN TEST1

author: Monllor-Berbegal, O.
"""

import numpy as np
from scipy.io import FortranFile as FF
import matplotlib.pyplot as plt
import sys
############################################################################################################
#IMPORT read_voids FROM masclet_framework: https://github.com/dvallesp/masclet_framework
############################################################################################################
sys.path.append('/home/monllor/projects')
from masclet_framework import read_voids
############################################################################################################

###############################
#open void mock catalogue
###############################

with open('mocks', 'r') as f:
    header = f.readline().split()
    nvoids = int(header[0])
    catalogue = []
    for i in range(nvoids):
        void = {}
        line = np.array(f.readline().split()).astype('f4')
        void['id'] = int(line[0])
        void['x'] = line[1]
        void['y'] = line[2]
        void['z'] = line[3]
        void['a'] = line[4]
        void['b'] = line[5]
        void['c'] = line[6]
        void['R'] = (void['a']*void['b']*void['c'])**(1./3.)
        catalogue.append(void)

    catalogue = {k: np.array([v[k] for v in catalogue]) for k in catalogue[0].keys()}

###############################
#now open voids found by AVISM
###############################


it = 1
path = '../output_files'
void_details, void_levels = read_voids.read_void_catalogue(it, path, 'arrays')
#void finder details
nlev = void_details['nlev']
levmin = void_details['levmin']
levmax = void_details['levmax']
nmax = void_details['nmax'] #levmin grid size !!!
nmay = void_details['nmay']
nmaz = void_details['nmaz']
L0 = void_details['L0']
#voids at each level
level = 0
ilevel = level-levmin
nvoids2 = void_levels[ilevel]['nvoids']
voids = void_levels[ilevel]['voids']
#void details
xc = voids['gxc'] #geometrical center
yc = voids['gyc']
zc = voids['gzc']
R = voids['R']
ellipticity = voids['ellipticity']
inv_porosity = voids['inv_porosity']

#keep only nvoids
nvoids2 = 1*nvoids
xc = xc[:nvoids2]
yc = yc[:nvoids2]
zc = zc[:nvoids2]
R = R[:nvoids2]
ellipticity = ellipticity[:nvoids2]
inv_porosity = inv_porosity[:nvoids2]


###############################
# match MOCK voids with voids found
###############################

avism2mock = np.zeros((nvoids,), dtype=int)
for ivoid in range(nvoids2):
    dx = xc[ivoid] - catalogue['x']
    dy = yc[ivoid] - catalogue['y']
    dz = zc[ivoid] - catalogue['z']
    dist = np.sqrt(dx**2 + dy**2 + dz**2)
    imatch = np.argmin(dist)
    avism2mock[ivoid] = imatch

###############################
#comparison MOCK voids vs voids found
###############################


#center offset
delta_CO = (   (xc - catalogue['x'][avism2mock])**2 
             + (yc - catalogue['y'][avism2mock])**2 
             + (zc - catalogue['z'][avism2mock])**2 )**0.5 / catalogue['R'][avism2mock]

#error in radius
delta_R = 1. - (R / catalogue['R'][avism2mock])

#error in ellipticity
delta_e = abs(ellipticity - catalogue['c'][avism2mock]/catalogue['a'][avism2mock])/catalogue['c'][avism2mock]

#error in inv_porosity
delta_ip = 1. - inv_porosity 

#divide in 3 bins: small, intermediate, large
mask_small = catalogue['R'][avism2mock] < 11.
mask_intermediate = (catalogue['R'][avism2mock] >= 11.) & (catalogue['R'][avism2mock] < 17.)
mask_large = catalogue['R'][avism2mock] >= 17.


#PLOT HISTOGRAMS OF ERRORS
fig, axs = plt.subplots(2,2, figsize=(7,4.5), dpi=200)

nbins = 20
mini = 0.
maxi = 1.2*np.max([np.max(delta_CO), np.max(delta_R), np.max(delta_e), np.max(delta_ip)])
width = 2
histtype = 'step'
bins = np.linspace(mini, maxi, nbins)
xlimlow = -0.002
xlimhigh = 0.352
ylimlow = 0.5
ylimhigh = 0.6*nvoids
alpha_mean = 0.8
width_mean = 1.
color_small = 'tab:cyan'
color_intermediate = 'tab:olive'
color_large = 'tab:red'

#center_offset
axs[0,0].hist(delta_CO[mask_small], bins=bins, alpha = 0.15, color=color_small)
axs[0,0].hist(delta_CO[mask_small], bins=bins, histtype=histtype, label='small', linewidth=width, color=color_small)
axs[0,0].hist(delta_CO[mask_intermediate], bins=bins, alpha = 0.15, color=color_intermediate)
axs[0,0].hist(delta_CO[mask_intermediate], bins=bins, histtype=histtype, label='intermediate', linewidth=width, color=color_intermediate)
axs[0,0].hist(delta_CO[mask_large], bins=bins, alpha = 0.15, color=color_large)
axs[0,0].hist(delta_CO[mask_large], bins=bins, histtype=histtype, label='large', linewidth=width, color=color_large)

axs[0,0].vlines(x=np.mean(delta_CO[mask_small]), color=color_small, linestyle='--', ymin=ylimlow, ymax=ylimhigh, alpha = alpha_mean, linewidth=width_mean)
axs[0,0].vlines(x=np.mean(delta_CO[mask_intermediate]), color=color_intermediate, linestyle='--', ymin=ylimlow, ymax=ylimhigh, alpha = alpha_mean, linewidth=width_mean)
axs[0,0].vlines(x=np.mean(delta_CO[mask_large]), color=color_large, linestyle='--', ymin=ylimlow, ymax=ylimhigh, alpha = alpha_mean, linewidth=width_mean)

axs[0,0].text(0.6, 0.4, r'$\langle \Delta_{CO} \rangle$' + f' = {np.mean(delta_CO[mask_small]):.3f}', transform=axs[0,0].transAxes, color=color_small, fontsize='small')
axs[0,0].text(0.6, 0.3, r'$\langle \Delta_{CO} \rangle$' + f' = {np.mean(delta_CO[mask_intermediate]):.3f}', transform=axs[0,0].transAxes, color=color_intermediate, fontsize='small')
axs[0,0].text(0.6, 0.2, r'$\langle \Delta_{CO} \rangle$' + f' = {np.mean(delta_CO[mask_large]):.3f}', transform=axs[0,0].transAxes, color=color_large, fontsize='small')

axs[0,0].set_xlabel('$\Delta_{CO}$')
axs[0,0].set_xlim(xlimlow, xlimhigh)
axs[0,0].set_ylim(ylimlow, ylimhigh)

#radius
axs[0,1].hist(delta_R[mask_small], bins=bins, alpha = 0.15, color=color_small)
axs[0,1].hist(delta_R[mask_small], bins=bins, histtype=histtype, label='small', linewidth=width, color=color_small)
axs[0,1].hist(delta_R[mask_intermediate], bins=bins, alpha = 0.15, color=color_intermediate)
axs[0,1].hist(delta_R[mask_intermediate], bins=bins, histtype=histtype, label='intermediate', linewidth=width, color=color_intermediate)
axs[0,1].hist(delta_R[mask_large], bins=bins, alpha = 0.15, color=color_large)
axs[0,1].hist(delta_R[mask_large], bins=bins, histtype=histtype, label='large', linewidth=width, color=color_large)

axs[0,1].vlines(x=np.mean(delta_R[mask_small]), color=color_small, linestyle='--', ymin=ylimlow, ymax=ylimhigh, alpha = alpha_mean, linewidth=width_mean)
axs[0,1].vlines(x=np.mean(delta_R[mask_intermediate]), color=color_intermediate, linestyle='--', ymin=ylimlow, ymax=ylimhigh, alpha = alpha_mean, linewidth=width_mean)
axs[0,1].vlines(x=np.mean(delta_R[mask_large]), color=color_large, linestyle='--', ymin=ylimlow, ymax=ylimhigh, alpha = alpha_mean, linewidth=width_mean)

axs[0,1].text(0.6, 0.4, r'$\langle \Delta_{R_e} \rangle$' + f' = {np.mean(delta_R[mask_small]):.3f}', transform=axs[0,1].transAxes, color=color_small, fontsize='small')
axs[0,1].text(0.6, 0.3, r'$\langle \Delta_{R_e} \rangle$' + f' = {np.mean(delta_R[mask_intermediate]):.3f}', transform=axs[0,1].transAxes, color=color_intermediate, fontsize='small')
axs[0,1].text(0.6, 0.2, r'$\langle \Delta_{R_e} \rangle$' + f' = {np.mean(delta_R[mask_large]):.3f}', transform=axs[0,1].transAxes, color=color_large, fontsize='small')

axs[0,1].set_xlabel('$\Delta_{R_e}$')
axs[0,1].set_xlim(xlimlow, xlimhigh)
axs[0,1].set_ylim(ylimlow, ylimhigh)

#ellipticity
axs[1,0].hist(delta_e[mask_small], bins=bins, alpha = 0.15, color=color_small)
axs[1,0].hist(delta_e[mask_small], bins=bins, histtype=histtype, label='small', linewidth=width, color=color_small)
axs[1,0].hist(delta_e[mask_intermediate], bins=bins, alpha = 0.15, color=color_intermediate)
axs[1,0].hist(delta_e[mask_intermediate], bins=bins, histtype=histtype, label='intermediate', linewidth=width, color=color_intermediate)
axs[1,0].hist(delta_e[mask_large], bins=bins, alpha = 0.15, color=color_large)
axs[1,0].hist(delta_e[mask_large], bins=bins, histtype=histtype, label='large', linewidth=width, color=color_large)

axs[1,0].vlines(x=np.mean(delta_e[mask_small]), color=color_small, linestyle='--', ymin=ylimlow, ymax=ylimhigh, alpha = alpha_mean, linewidth=width_mean)
axs[1,0].vlines(x=np.mean(delta_e[mask_intermediate]), color=color_intermediate, linestyle='--', ymin=ylimlow, ymax=ylimhigh, alpha = alpha_mean, linewidth=width_mean)
axs[1,0].vlines(x=np.mean(delta_e[mask_large]), color=color_large, linestyle='--', ymin=ylimlow, ymax=ylimhigh, alpha = alpha_mean, linewidth=width_mean)

axs[1,0].text(0.6, 0.4, r'$\langle \Delta_{\epsilon} \rangle$' + f' = {np.mean(delta_e[mask_small]):.3f}', transform=axs[1,0].transAxes, color=color_small, fontsize='small')
axs[1,0].text(0.6, 0.3, r'$\langle \Delta_{\epsilon} \rangle$' + f' = {np.mean(delta_e[mask_intermediate]):.3f}', transform=axs[1,0].transAxes, color=color_intermediate, fontsize='small')
axs[1,0].text(0.6, 0.2, r'$\langle \Delta_{\epsilon} \rangle$' + f' = {np.mean(delta_e[mask_large]):.3f}', transform=axs[1,0].transAxes, color=color_large, fontsize='small')

axs[1,0].set_xlabel('$\Delta_{\epsilon}$')
axs[1,0].set_xlim(xlimlow, xlimhigh)
axs[1,0].set_ylim(ylimlow, ylimhigh)

#inv_porosity
axs[1,1].hist(delta_ip[mask_small], bins=bins, alpha = 0.15, color=color_small)
axs[1,1].hist(delta_ip[mask_small], bins=bins, histtype=histtype, label='small', linewidth=width, color=color_small)
axs[1,1].hist(delta_ip[mask_intermediate], bins=bins, alpha = 0.15, color=color_intermediate)
axs[1,1].hist(delta_ip[mask_intermediate], bins=bins, histtype=histtype, label='intermediate', linewidth=width, color=color_intermediate)
axs[1,1].hist(delta_ip[mask_large], bins=bins, alpha = 0.15, color=color_large)
axs[1,1].hist(delta_ip[mask_large], bins=bins, histtype=histtype, label='large', linewidth=width, color=color_large)

axs[1,1].vlines(x=np.mean(delta_ip[mask_small]), color=color_small, linestyle='--', ymin=ylimlow, ymax=ylimhigh, alpha = alpha_mean, linewidth=width_mean)
axs[1,1].vlines(x=np.mean(delta_ip[mask_intermediate]), color=color_intermediate, linestyle='--', ymin=ylimlow, ymax=ylimhigh, alpha = alpha_mean, linewidth=width_mean)
axs[1,1].vlines(x=np.mean(delta_ip[mask_large]), color=color_large, linestyle='--', ymin=ylimlow, ymax=ylimhigh, alpha = alpha_mean, linewidth=width_mean)

axs[1,1].text(0.6, 0.4, r'$\langle \Delta_{IP} \rangle$' + f' = {np.mean(delta_ip[mask_small]):.3f}', transform=axs[1,1].transAxes, color=color_small, fontsize='small')
axs[1,1].text(0.6, 0.3, r'$\langle \Delta_{IP} \rangle$' + f' = {np.mean(delta_ip[mask_intermediate]):.3f}', transform=axs[1,1].transAxes, color=color_intermediate, fontsize='small')
axs[1,1].text(0.6, 0.2, r'$\langle \Delta_{IP} \rangle$' + f' = {np.mean(delta_ip[mask_large]):.3f}', transform=axs[1,1].transAxes, color=color_large, fontsize='small')

axs[1,1].set_xlabel('$\Delta_{IP}$')
axs[1,1].set_xlim(xlimlow, xlimhigh)
axs[1,1].set_ylim(ylimlow, ylimhigh)

#legend
fig.supylabel('Counts')
plt.legend(loc='upper right')

#y-scale
for ax in axs.ravel():
    ax.set_yscale('log')

#outward ticks
for ax in axs.ravel():
    ax.tick_params(which='both', direction='out', labelsize='large')
    ax.tick_params(which='major', length=4, width=1)
    ax.tick_params(which='minor', length=2, width=1)
    ax.set_xticks([0., 0.1, 0.2, 0.3])
    ax.minorticks_on()

for axis in ['top','bottom','left','right']:
    for ax in axs.ravel():
        ax.spines[axis].set_linewidth(1.3)

#remove gap
fig.tight_layout()

plt.savefig('test1_hist.png')



