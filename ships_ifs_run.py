############################
# Date: 09/04/2019
# Title: Running script for SHIPS for IFS data
# Description: Use this script to run SHIPS for IFS data. In this script you'll find all the necessary parameters to run SHIPS. ONLY SPHERE-DC DATA FOR NOW. VIP and pyKLIP are used.
# VIP version: 0.7.0 (Rainot edit.)
# pyKLIP version: 1.1
############################

# Set up your parameters

## Define images to analyse
cube_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/ifs_sortframes_dc-IFS_SCIENCE_REDUCED_SPECTRAL_MASTER_CUBE_SORTED-center_im_sorted.fits'
wavelength_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/ifs_sortframes_dc-IFS_SCIENCE_LAMBDA_INFO-lam.fits'
angles_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/ifs_sortframes_dc-IFS_SCIENCE_PARA_ROTATION_CUBE_SORTED-rotnth_sorted.fits'
psf_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/corrected_psf.fits'

## Setup Parameters
fwhm = 4. # FWHM default value: 4.

## Photometry
comp_pos = (112,55) # Companion position in pixels (X,Y)
psf_pos = (32, 33) # PSF position in pixels (X,Y)
radial_dist = 97.684904197311155 # Radial distance of companion in pixels
position_angle = 159.37210826761003  # Position angle of companion in degrees
noise_aperture_pos_comp = (92,102) # Position in pixels of the circular annulus aperture for noise measurement in the case of the companion
noise_aperture_pos_psf = (12,22) # Position in pixels of the circular annulus aperture for noise measurement in the case of the PSF

## PCA
ncomp_pca = 1 # Number of principal components for PCA

## Do you want to see the image?
see_cube = True
see_collapsed_cube = True



## Contrast curves
contrast_curves = True # True or False
n_branches = 1 # Number of branches for contrast curves


# ---------------------------------------------------------------------------

# Running script (DO NOT MODIFY)

# Some definitions

## Load libraries
import __init__
from matplotlib.pyplot import *
import vip
from vip.phot.fakecomp import inject_fcs_cube, inject_fc_frame, psf_norm
import numpy as np
import scipy
import astropy.io.fits as fits
from astropy import wcs
import photutils
from photutils import CircularAperture
from photutils import CircularAnnulus
import matplotlib.pyplot as plt
import glob
import pyklip.instruments.SPHERE as SPHERE
import math as mh
import pyklip.parallelized as parallelized
from scipy.integrate import quad, dblquad

## Define constants
c = 299792458. # Speed of light
Ro = 6.957e8 # Solar Radius
sr2pc = 44334448.0068964 # Convert steraradians to parsec
pxscale = 0.0074 # IFS pixel scale in arcsec/pixel
PA = position_angle + 90 # Correct for VIP unconventional rotation axis
rad_dists = [radial_dist] # Convert radial distance float to array, otherwise negfc won't work

## Open image files
cube = vip.fits.open_fits(cube_filepath)
wl = vip.fits.open_fits(wavelength_filepath)
angs = vip.fits.open_fits(angles_filepath)
psf = vip.fits.open_fits(psf_filepath)

## Define some Parameters
psf_scaled = np.zeros_like(psf) # The psf will need to be scaled
flevel = np.zeros_like(cube[:,0,0,0]) # Flux level for the companion
flevel = np.array(flevel) # Redefinition - why?

# Stellar photometry of the companion

## Collapse the images for better photometry measurement
cube_derot = vip.preproc.cube_derotate(cube,angs) # Rotate the images to the same north
cube_wl_coll = vip.preproc.cube_collapse(cube_derot,wl_cube=True) # Collapse along the rotation axis - 3D image
cube_coll = vip.preproc.cube_collapse(cube_derot,wl_cube=False) # Collapse along the wavelength axis - 2D image

## Aperture photometry of companions and PSF

### Define photometry
noise_phot = np.zeros_like(wl) #Noise photometry
psf_final_sum = np.zeros_like(wl) #PSF photometry
final_sum = np.zeros_like(wl) #Companion photometry

### Apertures
aper_comp = photutils.CircularAperture(comp_pos, 5./2*fwhm)
aper_psf = photutils.CircularAperture(psf_pos, 5./2*fwhm)
aper_noise_comp = photutils.CircularAnnulus((145,145),noise_aperture_pos_comp[0],noise_aperture_pos_comp[1])
aper_noise_psf = photutils.CircularAnnulus(psf_pos,noise_aperture_pos_psf[0],noise_aperture_pos_psf[1])

### Aperture photometry
for i in range(0,wl.shape[0]):
    #Noise
    phot_noise = photutils.aperture_photometry(cube_wl_coll[i], aper_noise_comp)
    noise_phot[i] = np.array(phot_noise['aperture_sum'])
    #PSF
    phot_psf = photutils.aperture_photometry(psf[i], aper_psf)
    phot_psf_noise = photutils.aperture_photometry(psf[i], aper_noise_psf)
    psf_bkg_mean = phot_psf_noise['aperture_sum'] / aper_noise_psf.area()
    psf_bkg_sum = psf_bkg_mean * aper_psf.area()
    psf_final_sum[i] = phot_psf['aperture_sum'] - psf_bkg_sum
    #Companion
    phot = photutils.aperture_photometry(cube_wl_coll[i], aper_comp)
    bkg_mean = (phot_noise['aperture_sum']-phot['aperture_sum']) / (aper_noise_comp.area()-aper_comp.area())
    bkg_sum = bkg_mean * aper_comp.area()
    final_sum[i] = phot['aperture_sum'] - bkg_sum

## Contrast curve
if contrast_curves == True:
    #Crop the PSF to match the size of companion
    for i in range(psf.shape[0]):
                psf_temp = vip.phot.psf_norm(psf[i], size=3*fwhm, fwhm=fwhm)
                if i==0:
                    psf_crop = np.zeros((len(psf),psf_temp.shape[0],psf_temp.shape[1]))
                psf_crop[i] = psf_temp
    cube_negfc = vip.phot.inject_fcs_cube(cube,psf_crop,-angs,flevel=-105,plsc=pxscale,rad_dists=rad_dists,theta=PA) # Remove companion using NEGFC technique

    print("Companion removed")
    print("Computing contrast curve...")
    contrcurve = vip.phot.contrast_curve(cube_negfc,-angs,psf,4.,pxscale,psf_final_sum,vip.pca.pca,nbranch=n_branches,
               dpi=300, student=False, debug=True ,plot=True, verbose=True, ncomp=ncomp_pca, scale_list=wl)

elif contrast_curves == False:
    print("No contrast curve")
