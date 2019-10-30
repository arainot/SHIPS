############################
# Date: 09/10/2019
# Title: Auto running script for SHIPS for IFS data for the MCMC computation. Requires the simplex filefrom SHIPS.
# Description: Use this script to run the MCMC computation of SHIPS in auto mode for IFS data. Requires the simplex filefrom SHIPS. ONLY SPHERE-DC DATA FOR NOW. VIP is used.
# VIP version: 0.9.11
# Python version: 3 ONLY
############################

# Set up your parameters

## Read thewavelength number
import sys
i = sys.argv[1] # read the wavelength index number
i = int(i)

## Define images to analyse
wavelength_filepath = '/home/alan/data/Backup_macbook/SPHERE/IFS/HD93403/ifs_sortframes_dc-IFS_SCIENCE_LAMBDA_INFO-lam.fits'
cube_filepath = '/home/alan/data/Backup_macbook/SPHERE/IFS/HD93403/ifs_sortframes_dc-IFS_SCIENCE_REDUCED_SPECTRAL_MASTER_CUBE_SORTED-center_im_sorted.fits'
angles_filepath = '/home/alan/data/Backup_macbook/SPHERE/IFS/HD93403/ifs_sortframes_dc-IFS_SCIENCE_PARA_ROTATION_CUBE_SORTED-rotnth_sorted.fits'
psf_filepath = '/home/alan/data/Backup_macbook/SPHERE/IFS/HD93403/ifs_sortframes_dc-IFS_SCIENCE_PSF_MASTER_CUBE-median_unsat.fits'

## Photometry
comp_pos = (112.,54.) # Companion position in pixels from the center of the frame (X,Y)
psf_pos = (32, 33) # PSF position in pixels (X,Y)
radial_dist = 97. # Radial distance of companion in pixels
position_angle = 159.  # Position angle of companion in degrees
noise_aperture_pos_comp = (92,102) # Position in pixels of the circular annulus aperture for noise measurement in the case of the companion
noise_aperture_pos_psf = (12,22) # Position in pixels of the circular annulus aperture for noise measurement in the case of the PSF

## Computing power
ncores = 4 # Number of cores you are willing to share for the computation

## PCA
ncomp_pca = 1 # Number of principal components for PCA

## Spectrum extraction with Simplex Nelder-Mead optimisation
sspec_file = '/home/alan/data/Backup_macbook/SPHERE/IFS/HD93403/VIP_simplex.txt' # Filepath to save the Simplex spectrum

## Spectrum extraction with MCMC
extract_mcmc = True # Will compute the MCMC for all 39 wavelengths !! This takes ~1,5h per wavelength and is very computer intensive !!
ann_width = 3 # Annulus width of Simplex
aper_radius = 3 # Aperture Radius of PCA
source = 'HD93403' # Give name for your source
mcmc_path = '/home/alan/Documents/Thesis/SPHERE/spectra/HD93403/mcmc11/' # Directory where MCMC results will be stored
plot_mcmc = False # Plot the mcmc errors with simplex?

# ---------------------------------------------------------------------------

# Running script (DO NOT MODIFY)

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
import photutils
from photutils import CircularAperture
from photutils import CircularAnnulus
import matplotlib.pyplot as plt
import glob
#import pyklip.instruments.SPHERE as SPHERE
import math as mh
#import pyklip.parallelized as parallelized
from scipy.integrate import quad, dblquad

## Define constants
c = 299792458. # Speed of light
Ro = 6.957e8 # Solar Radius
sr2pc = 44334448.0068964 # Convert steraradians to parsec
pxscale = 0.0074 # IFS pixel scale in arcsec/pixel
PA = position_angle + 90 # Correct for VIP unconventional rotation axis

## Open image files
cube = vip_hci.fits.open_fits(cube_filepath)
wl = vip_hci.fits.open_fits(wavelength_filepath)
angs = vip_hci.fits.open_fits(angles_filepath)
psf = vip_hci.fits.open_fits(psf_filepath)

## Define some Parameters
psf = np.median(psf, axis=1) # Take the median of all PSFs
psf_scaled = np.zeros_like(psf) # The psf will need to be scaled
flevel = np.zeros_like(cube[:,0,0,0]) # Flux level for the companion
flevel = np.array(flevel) # Redefinition - why?

## Get FWHM of images & normalised PSF
psf_norm, maxflux, fwhm = vip_hci.metrics.normalize_psf(psf, fwhm='fit', size=31,verbose=False,full_output=True) # maxflux is a dummy variable

# # Stellar photometry of the companion
#
# ## Collapse the images for better photometry measurement
# cube_derot = vip_hci.preproc.cube_derotate(cube,angs) # Rotate the images to the same north
# cube_wl_coll = vip_hci.preproc.cube_collapse(cube_derot,wl_cube=True) # Collapse along the rotation axis - 3D image
# cube_coll = vip_hci.preproc.cube_collapse(cube_derot,wl_cube=False) # Collapse along the wavelength axis - 2D image
#
# ## Aperture photometry of companions and PSF
#
# ### Define photometry
# noise_phot = np.zeros_like(wl) #Noise photometry
# psf_final_sum = np.zeros_like(wl) #PSF photometry
# final_sum = np.zeros_like(wl) #Companion photometry
#
# ### Apertures
# aper_noise_comp = photutils.CircularAnnulus((145,145),noise_aperture_pos_comp[0],noise_aperture_pos_comp[1])
# aper_noise_psf = photutils.CircularAnnulus(psf_pos,noise_aperture_pos_psf[0],noise_aperture_pos_psf[1])
#
# ### Aperture photometry
# for i in range(0,wl.shape[0]):
#     ### Apertures dependent on channel
#     aper_comp = photutils.CircularAperture(comp_pos, 1./2*fwhm[i])
#     aper_psf = photutils.CircularAperture(psf_pos, 1./2*fwhm[i])
#     ### Noise
#     phot_noise = photutils.aperture_photometry(cube_wl_coll[i], aper_noise_comp)
#     noise_phot[i] = np.array(phot_noise['aperture_sum'])
#     ### PSF
#     phot_psf = photutils.aperture_photometry(psf[i], aper_psf)
#     phot_psf_noise = photutils.aperture_photometry(psf[i], aper_noise_psf)
#     psf_bkg_mean = phot_psf_noise['aperture_sum'] / aper_noise_psf.area
#     psf_bkg_sum = psf_bkg_mean * aper_psf.area
#     psf_final_sum[i] = phot_psf['aperture_sum'] - psf_bkg_sum
#     ### Companion
#     phot = photutils.aperture_photometry(cube_wl_coll[i], aper_comp)
#     bkg_mean = (phot_noise['aperture_sum']-phot['aperture_sum']) / (aper_noise_comp.area-aper_comp.area)
#     bkg_sum = bkg_mean * aper_comp.area
#     final_sum[i] = phot['aperture_sum'] - bkg_sum

# Spectrum extraction with NM
f = open(sspec_file,'r') # Read the spectrum file
simplex_guess = np.zeros((39,3)) # Set the simplex variable: r, PA, flux
j=0 # Iterator for the simplex variable
for line in f:
    line = line.strip()
    columns = line.split()
    simplex_guess[j][0] = float(columns[0])
    simplex_guess[j][1] = float(columns[1])
    simplex_guess[j][2] = float(columns[2])
    j+=1
f.close()
print("Spectrum successfully loaded!")

# Spectrum extraction with MCMC
if extract_mcmc == True:
    instru= 'IFS36059' # Define instrument parameters
    fig_merit='sum' # Summation figure of merit
    outpath = mcmc_path.format(source) # Path to save MCMC files

    print("########## MCMC Sampling starting... ##########")

    nwalkers, itermin, itermax = (100,200,500) # as recommended by Oli W

    initialState = simplex_guess[i] # Take r, PA and flux from simplex
    bounds=[[0.75*initialState[0],1.25*initialState[0]],[0.75*initialState[1],1.25*initialState[1]],[0.75*initialState[2],1.30*initialState[2]]] # Initiate bounds
    output_file = source+'_IFS_wavelength_{}'.format(i) # Save to output file

    chain_40 = vip_hci.negfc.mcmc_negfc_sampling(cube[i], -angs,  psf_norm[i], ncomp_pca, pxscale, initialState, ann_width,
                                             aper_radius, cube_ref=None, svd_mode='lapack', nwalkers=nwalkers,
                                             bounds=bounds, niteration_min=itermin,
                                             niteration_limit=itermax, check_maxgap=50, nproc= ncores,
                                             output_file=output_file, display=False,verbosity=1, save=True,
                                             rhat_threshold=1.01, niteration_supp=0, fmerit=fig_merit) # MCMC run per channel
    print("########## MCMC Sampling done! ##########")
