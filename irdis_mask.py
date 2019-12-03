############################
# Date: 26/11/2019
# Title: Annulus mask
# Description: Use this script to mask the data in IRDIS by calculating the noise in an annular fashion and replacing the companion by that noise mask.
# VIP version: 0.9.11
# Python version: 3
############################

# Set up your parameters

## Define images to analyse
cube_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/ird_convert_dc-IRD_SCIENCE_REDUCED_MASTER_CUBE-center_im.fits'
cube_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/cube_negfc_new.fits'
wavelength_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/ird_convert_dc-IRD_SCIENCE_LAMBDA_INFO-lam.fits'
angles_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/ird_convert_dc-IRD_SCIENCE_PARA_ROTATION_CUBE-rotnth.fits'
psf_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/ird_convert_recenter_dc5-IRD_SCIENCE_PSF_MASTER_CUBE-median_unsat.fits'
# wavelength_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/Hugues_data/IRDIS/CEN3/ird_convert_recenter_dc-IRD_SCIENCE_LAMBDA_INFO-lam.fits'
# cube_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/Hugues_data/IRDIS/CEN3/ird_convert_recenter_dc-IRD_SCIENCE_REDUCED_MASTER_CUBE-center_im.fits'
# angles_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/Hugues_data/IRDIS/CEN3/ird_convert_recenter_dc-IRD_SCIENCE_PARA_ROTATION_CUBE-rotnth.fits'
# psf_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/Hugues_data/IRDIS/CEN3/ird_convert_recenter_dc-IRD_SCIENCE_PSF_MASTER_CUBE-median_unsat.fits'

## Photometry
comp_pos = ([491,456],[559,579],[688,630],[570,702],[311,486],[696,419],[296,472],[654,707],[519,243],[227,569],[346,778],[847,491],[899,507],[72,533],[819,180],[451,60],[49,396],[732,44],[648,16]) # Companion position in pixels (X,Y)
psf_pos = (33, 33) # PSF position in pixels (X,Y)
radial_dist = [60.16643583,81.60882305,210.78899402,197.121377136,201.68291946,211,219.,240.63665556,269.09106265,290.11,313.16,334.845228417,384.760748992,442.00180288,451.91383567,456.697517309,477.277760379,515.37656136072,515.7130985344468] # Radial distance of companion in pixels
position_angle = [295.9,311.8] # Position angle of companion in degrees
noise_aperture_pos_comp = (512,512) # Position in pixels of the circular annulus aperture for noise measurement in the case of the companion
noise_aperture_pos_psf = (12,22) # Position in pixels of the circular annulus aperture for noise measurement in the case of the PSF

## Computing power
ncores = 4 # Number of cores you are willing to share for the computation

## Do you want to see the image?
see_cube = False # Original cube
see_collapsed_cube = False # Collapsed cube
see_psf_norm = False # Normalised PSF
see_cube_centre = False # Check if the image is centered correctly
size_psf = 31


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
import random
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

# Stellar photometry of the companion

## Collapse the images for better photometry measurement
cube_wl_coll = np.zeros_like(cube[:,0,:,:])
for i in range(len(wl)):
        cube_wl_coll[i] = vip_hci.hci_postproc.median_sub(cube[i],-angs,fwhm=fwhm[i],verbose=False) # Rotate & collapse along the rotation axis - 3D image
cube_derot = vip_hci.preproc.cube_derotate(cube,angs) # Rotate the images to the same north
# cube_wl_coll = vip_hci.preproc.cube_collapse(cube_derot,wl_cube=True) # Collapse along the rotation axis - 3D image
#cube_coll = vip_hci.preproc.cube_collapse(cube_derot,wl_cube=False) # Collapse along the wavelength axis - 2D image

## Aperture photometry of companions and PSF

### Define photometry
noise_phot = np.zeros_like(wl) #Noise photometry
psf_final_sum = np.zeros_like(wl) #PSF photometry
final_sum_K1 = np.zeros_like(radial_dist) #Companion photometry in the K1 band
final_sum_K2 = np.zeros_like(radial_dist) #Companion photometry in the K2 band

### Apertures
aper_noise_psf = photutils.CircularAnnulus(psf_pos,noise_aperture_pos_psf[0],noise_aperture_pos_psf[1])

### Aperture photometry - PSF
for i in range(0,len(wl)):
    ### Aperture
    aper_psf = photutils.CircularAperture(psf_pos, 1./2*fwhm[i])
    ### Flux
    phot_psf = photutils.aperture_photometry(psf[i], aper_psf)
    phot_psf_noise = photutils.aperture_photometry(psf[i], aper_noise_psf)
    psf_bkg_mean = phot_psf_noise['aperture_sum'] / aper_noise_psf.area()
    psf_bkg_sum = psf_bkg_mean * aper_psf.area()
    psf_final_sum[i] = phot_psf['aperture_sum'] - psf_bkg_sum

### Aperture photometry - Companions
for i in range(0,len(radial_dist)):
    ### Apertures dependent on companions
    aper_noise_comp = photutils.CircularAnnulus((512,512),radial_dist[i]-5,radial_dist[i]+5)
    aper_comp_K1 = photutils.CircularAperture((comp_pos[i][0], comp_pos[i][1]),1./2*fwhm[0])
    aper_comp_K2 = photutils.CircularAperture((comp_pos[i][0], comp_pos[i][1]),1./2*fwhm[1])
    ### Flux
    phot_noise_K1 = photutils.aperture_photometry(cube_wl_coll[0], aper_noise_comp)
    phot_noise_K2 = photutils.aperture_photometry(cube_wl_coll[1], aper_noise_comp)
    phot_K1 = photutils.aperture_photometry(cube_wl_coll[0], aper_comp_K1)
    phot_K2 = photutils.aperture_photometry(cube_wl_coll[1], aper_comp_K2)
    bkg_mean_K1 = (phot_noise_K1['aperture_sum']-phot_K1['aperture_sum']) / (aper_noise_comp.area()-aper_comp_K1.area())
    bkg_mean_K2 = (phot_noise_K2['aperture_sum']-phot_K2['aperture_sum']) / (aper_noise_comp.area()-aper_comp_K2.area())
    bkg_sum_K1 = bkg_mean_K1 * aper_comp_K1.area()
    bkg_sum_K2 = bkg_mean_K2 * aper_comp_K2.area()
    final_sum_K1[i] = phot_K1['aperture_sum'] - bkg_sum_K1
    final_sum_K2[i] = phot_K2['aperture_sum'] - bkg_sum_K2

cube_mask = cube # Make a safe copy of the cube that will be used as the masked cube

# Masking
## Ad
for i in range(512,542): # Iterate through all pixels in X/Y positions
    for j in range(560,587):
        X = random.randint(548,558) # Find random coordinates within a predefined box for noise search
        Y = random.randint(447,464)
        mask_value = cube[:,:,Y,X] # Find the noise values at the random coordinates
        cube_mask[:,:,j,i] = mask_value

## Ab
for i in range(435,506): # Iterate through all pixels in X/Y positions
    for j in range(409,467):
        X = random.randint(428,446) # Find random coordinates within a predefined box for noise search
        Y = random.randint(536,545)
        mask_value = cube[:,:,Y,X] # Find the noise values at the random coordinates
        cube_mask[:,:,j,i] = mask_value

## E
for i in range(339,377): # Iterate through all pixels in X/Y positions
    for j in range(351,392):
        X = random.randint(714,727) # Find random coordinates within a predefined box for noise search
        Y = random.randint(500,512)
        mask_value = cube[:,:,Y,X] # Find the noise values at the random coordinates
        cube_mask[:,:,j,i] = mask_value

## S15
for i in range(435,506): # Iterate through all pixels in X/Y positions
    for j in range(409,467):
        X = random.randint(428,446) # Find random coordinates within a predefined box for noise search
        Y = random.randint(536,545)
        mask_value = cube[:,:,Y,X] # Find the noise values at the random coordinates
        cube_mask[:,:,j,i] = mask_value

## S16
for i in range(339,377): # Iterate through all pixels in X/Y positions
    for j in range(351,392):
        X = random.randint(714,727) # Find random coordinates within a predefined box for noise search
        Y = random.randint(500,512)
        mask_value = cube[:,:,Y,X] # Find the noise values at the random coordinates
        cube_mask[:,:,j,i] = mask_value

d = vip_hci.hci_postproc.median_sub(cube_mask[0],-angs,fwhm=fwhm[0],verbose=False)
ds9 = vip_hci.Ds9Window()
ds9.display(cube_mask[0,0],d)
