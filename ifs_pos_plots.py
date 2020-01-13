############################
# Date: 20/12/2019
# Title: Script for plotting the positions of IFS sources
# Description: Use this script make a plot comparing the different values of position and errors from different sources
# VIP version: 0.9.11 (Rainot edit.)
# Python version: 3 ONLY
############################
# Some definitions

## Load libraries
import __init__
import sys
import matplotlib
import vip_hci
from hciplot import plot_frames, plot_cubes
from vip_hci.metrics.fakecomp import cube_inject_companions, frame_inject_companion, normalize_psf
import numpy as np
import scipy
import astropy.io.fits as fits
from astropy import wcs
import matplotlib.pyplot as plt
import math as mh

# Find the files containing the values
mcmc_pos_file = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/IFS_MCMC_pos.txt'
mcmc_PA = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/IFS_MCMC_PA.txt'
mcmc_flux = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/IFS_MCMC_flux.txt'

mcmc_pos = np.array([])
mcmc_PA = np.array([])
mcmc_flux = np.array([])
r_K2 = np.array([])
theta_K2 = np.array([])
f_K2 = np.array([])
std_r_K1 = np.zeros_like(final_sum_K1)
std_theta_K1 = np.zeros_like(final_sum_K1)
std_flux_K1 = np.zeros_like(final_sum_K1)
std_r_K2 = np.zeros_like(final_sum_K1)
std_theta_K2 = np.zeros_like(final_sum_K1)
std_flux_K2 = np.zeros_like(final_sum_K1)

with open(mcmc_pos_file) as f:
    header1 = f.readline()
    header2 = f.readline()
    header3 = f.readline()
    for line in f:
        line = line.strip()
        columns = line.split()
        r_K2 = np.append(r_K2,float(columns[0]))
        theta_K2 = np.append(theta_K2,float(columns[1]))
        f_K2 = np.append(f_K2,float(columns[2]))
f.close


plt.grid(True)
plt.legend()
plt.figure(figsize=(12, 9))
fig = plt.figure(figsize=(12, 9))
ax1 = fig.add_subplot(111)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(prop={'size': 20})
# Ensure that the axis ticks only show up on the bottom and left of the plot.
# Ticks on the right and top of the plot are generally unnecessary chartjunk.
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.gca().invert_yaxis()
# Limit the range of the plot to only where the data is.
# Avoid unnecessary whitespace.
plt.ylim(12, 6)
plt.xlim(0.1,0.85)
# Make sure your axis ticks are large enough to be easily read.
# You don't want your viewers squinting to read your plot.
plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.5)


ax1.set_xlabel('Angular separation ["]', fontsize=20)
ax1.set_ylabel("$\Delta$$H$", fontsize=20)
ax1.plot(contrast_c_dist[16:110], contrast_c_mags[16:110],lw=2.8,label='YJH',c='blue')
#ax1.plot(contrast_c_theory_dist[5:15]/1e3, contrast_c_theory_mags[5:15],lw=2.8,label='ETC',c='orange')
ax1.tick_params(axis='y')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set_ylim(ax1.get_ylim())
ax2.set_ylabel('Mass [$M_{\odot}$]',fontsize=18)  # we already handled the x-label with ax1
ax2.plot(contrast_c_dist[17:110]/1e3, contrast_c_mags[17:110], color="w",alpha=0)
ax2.set_yticks([6,7,8,9,10,11,12])
ax2.set_yticklabels(['0.1','0.4','0.7','1.4','2.0','4.0','7.0'],fontsize=18)
plt.xlim(0.1,0.9)
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()
plt.gca().invert_yaxis()
ax1.legend(prop={'size': 20})
plt.show()
#fig.savefig("/Users/alan/Nextcloud/PhD/Thesis/SPHERE/P96_CHIPS/Plots/Contrast_curves/QZCar/contrast_IFS.png",dpi=300,format='png')
