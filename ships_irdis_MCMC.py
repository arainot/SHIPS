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
cube_filepath = '/home/alan/data/Backup_macbook/SPHERE/IRDIS/QZCar/ird_convert_dc-IRD_SCIENCE_REDUCED_MASTER_CUBE-center_im.fits'
wavelength_filepath = '/home/alan/data/Backup_macbook/SPHERE/IRDIS/QZCar/ird_convert_dc-IRD_SCIENCE_LAMBDA_INFO-lam.fits'
angles_filepath = '/home/alan/data/Backup_macbook/SPHERE/IRDIS/QZCar/ird_convert_dc-IRD_SCIENCE_PARA_ROTATION_CUBE-rotnth.fits'
psf_filepath = '/home/alan/data/Backup_macbook/SPHERE/IRDIS/QZCar/ird_convert_recenter_dc5-IRD_SCIENCE_PSF_MASTER_CUBE-median_unsat.fits'

## Photometry
comp_pos = ([491,456],[559,579],[688,630],[570,702],[311,486],[696,419],[296,472],[654,707],[519,243],[227,569],[346,778],[847,491],[899,507],[72,533],[819,180],[451,60],[49,396],[732,44],[648,16]) # Companion position in pixels (X,Y)
psf_pos = (33, 33) # PSF position in pixels (X,Y)
radial_dist = [60.16643583,81.60882305,210.78899402,197.121377136,201.68291946,211,219.,240.63665556,269.09106265,290.11,313.16,334.845228417,384.760748992,442.00180288,451.91383567,456.697517309,477.277760379,515.37656136072,515.7130985344468] # Radial distance of companion in pixels
position_angle = [295.9,311.8] # Position angle of companion in degrees
noise_aperture_pos_comp = (512,512) # Position in pixels of the circular annulus aperture for noise measurement in the case of the companion
noise_aperture_pos_psf = (12,22) # Position in pixels of the circular annulus aperture for noise measurement in the case of the PSF

simplex_guess_K1 = ([334.845228417,266.28523872,11.6100343043],[384.760748992,269.212066163,11.6429408139])
simplex_guess_K2 = ([334.845228417,266.28523872,11.6100343043],[384.760748992,269.212066163,11.6429408139])

## Computing power
ncores = 4 # Number of cores you are willing to share for the computation

size_psf = 31

## PCA
ncomp_pca = 1 # Number of principal components for PCA

## Spectrum extraction with Simplex Nelder-Mead optimisation
ann_width = 3 # Annulus width of Simplex
aper_radius = 3 # Aperture Radius of PCA
#sspec_file_K1 = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/VIP_simplex_K1.txt' # Filepath to save the Simplex spectrum for the K1 band
#sspec_file_K2 = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/VIP_simplex_K2.txt' # Filepath to save the Simplex spectrum for the K2 band

## Reading MCMC results
extract_mcmc = True # Will compute the MCMC for all 39 wavelengths !! This takes ~1,5h per wavelength and is very computer intensive !!
ann_width = 3 # Annulus width of Simplex
aper_radius = 3 # Aperture Radius of PCA
source = 'QZCar' # Give name for your source
mcmc_path = '/home/alan/Desktop/' # Directory where MCMC results are stored


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
pxscale = 0.1225 # IRDIS pixel scale in arcsec/pixel
PA = np.array(position_angle) + 90 # Correct for VIP unconventional rotation axis

## Open image files
cube = vip_hci.fits.open_fits(cube_filepath)
wl = vip_hci.fits.open_fits(wavelength_filepath)
angs = vip_hci.fits.open_fits(angles_filepath)
psf = vip_hci.fits.open_fits(psf_filepath)

## Define some Parameters
psf = np.median(psf, axis=1) # Take the median value of all psf observations
psf_scaled = np.zeros_like(psf) # The psf will need to be scaled
flevel = np.zeros_like(cube[:,0,0,0]) # Flux level for the companion
flevel = np.array(flevel) # Redefinition - why?

## Get FWHM of images & normalised PSF
psf_med = vip_hci.preproc.cosmetics.cube_crop_frames(psf, size_psf, xy=(32, 32), verbose=True, force=True) # Resize the PSF
psf_norm, maxflux, fwhm = vip_hci.metrics.normalize_psf(psf_med, fwhm='fit',size=None, threshold=None,mask_core=None, model='gauss',imlib='opencv',interpolation='lanczos4',force_odd=True,full_output=True,verbose=False) # maxflux is a dummy variable

# Spectrum extraction with NM
# ##K1
# f = open(sspec_file_K1,'r') # Read the spectrum file
# simplex_guess_K1 = np.zeros((39,3)) # Set the simplex variable: r, PA, flux
# j=0 # Iterator for the simplex variable
# for line in f:
#     line = line.strip()
#     columns = line.split()
#     simplex_guess_K1[j][0] = float(columns[0])
#     simplex_guess_K1[j][1] = float(columns[1])
#     simplex_guess_K1[j][2] = float(columns[2])
#     j+=1
# f.close()
#
# ##K2
# f = open(sspec_file_K2,'r') # Read the spectrum file
# simplex_guess_K2 = np.zeros((39,3)) # Set the simplex variable: r, PA, flux
# j=0 # Iterator for the simplex variable
# for line in f:
#     line = line.strip()
#     columns = line.split()
#     simplex_guess_K2[j][0] = float(columns[0])
#     simplex_guess_K2[j][1] = float(columns[1])
#     simplex_guess_K2[j][2] = float(columns[2])
#     j+=1
# f.close()
# print("Spectrum successfully loaded!")

# Spectrum extraction with MCMC
if extract_mcmc == True:
    instru= 'IRDIS36059' # Define instrument parameters
    ann_width=annulus_width # Annulus width of MCMC
    aperture_radius=aperture_width # Aperture radius
    fig_merit='sum' # Summation figure of merit
    outpath = mcmc_path.format(source) # Path to save MCMC files

    print("########## MCMC Sampling starting... ##########")

    nwalkers, itermin, itermax = (150,200,1200) # as recommended by Oli W
    #K1
    initialState = simplex_guess_K1[i] # Take r, PA and flux from simplex
    bounds=[[0.75*initialState[0],1.25*initialState[0]],[0.75*initialState[1],1.25*initialState[1]],[0.75*initialState[2],1.30*initialState[2]]] # Initiate bounds
    output_file = source+'K1_IRDIS_companions_S{}'.format(i) # Save to output file

    chain_40 = vip.negfc.mcmc_negfc_sampling(cube[0], -angs,  psf_norm[0], ncomp_pca, pxscale, initialState, ann_width,
                                             aperture_radius, cube_ref=None, svd_mode='lapack', nwalkers=nwalkers,
                                             bounds=bounds, niteration_min=itermin,
                                             niteration_limit=itermax, check_maxgap=50, nproc= ncores,
                                             output_file=output_file, display=True,verbose=True, save=True,
                                             rhat_threshold=1.01, niteration_supp=0, fmerit=fig_merit) # MCMC run per channel
    #K2
    for i in range(0,2): # For each companion
    initialState = simplex_guess_K2[i] # Take r, PA and flux from simplex
    bounds=[[0.75*initialState[0],1.25*initialState[0]],[0.75*initialState[1],1.25*initialState[1]],[0.75*initialState[2],1.30*initialState[2]]] # Initiate bounds
    output_file = source+'K2_IRDIS_companions_S{}'.format(i) # Save to output file

    chain_40 = vip.negfc.mcmc_negfc_sampling(cube[1], -angs,  psf_scaled[1], ncomp_pca, pxscale, initialState, ann_width,
                                             aperture_radius, cube_ref=None, svd_mode='lapack', nwalkers=nwalkers,
                                             bounds=bounds, niteration_min=itermin,
                                             niteration_limit=itermax, check_maxgap=50, nproc= ncores,
                                             output_file=output_file, display=True,verbose=True, save=True,
                                             rhat_threshold=1.01, niteration_supp=0, fmerit=fig_merit) # MCMC run per channel
    print("########## MCMC Sampling done! ##########")
