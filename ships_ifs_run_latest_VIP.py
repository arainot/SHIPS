############################
# Date: 09/04/2019
# Title: Running script for SHIPS for IFS data
# Description: Use this script to run SHIPS for IFS data. In this script you'll find all the necessary parameters to run SHIPS. ONLY SPHERE-DC DATA FOR NOW. VIP and pyKLIP are used.
# VIP version: 0.9.8 (Rainot edit.)
# pyKLIP version: 1.1 NOT IMPLEMENTED YET
############################

# Set up your parameters

## Define images to analyse
cube_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/ifs_sortframes_dc-IFS_SCIENCE_REDUCED_SPECTRAL_MASTER_CUBE_SORTED-center_im_sorted.fits'
wavelength_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/ifs_sortframes_dc-IFS_SCIENCE_LAMBDA_INFO-lam.fits'
angles_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/ifs_sortframes_dc-IFS_SCIENCE_PARA_ROTATION_CUBE_SORTED-rotnth_sorted.fits'
psf_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/corrected_psf.fits'

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
see_cube = False # Original cube
see_collapsed_cube = False # Collapsed cube

## Contrast curves
contrast_curves = False # True or False
n_branches = 1 # Number of branches for contrast curves


# ---------------------------------------------------------------------------

# Running script (DO NOT MODIFY)

# Some definitions

## Load libraries
import __init__
import matplotlib
matplotlib.use('PS')
from matplotlib.pyplot import *
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
rad_dists = [radial_dist] # Convert radial distance float to array, otherwise negfc won't work

## Open image files
cube = vip_hci.fits.open_fits(cube_filepath)
wl = vip_hci.fits.open_fits(wavelength_filepath)
angs = vip_hci.fits.open_fits(angles_filepath)
psf = vip_hci.fits.open_fits(psf_filepath)

## Define some Parameters
psf_scaled = np.zeros_like(psf) # The psf will need to be scaled
flevel = np.zeros_like(cube[:,0,0,0]) # Flux level for the companion
flevel = np.array(flevel) # Redefinition - why?

## Check the RAW data cubes
if see_cube == True:
    ds9 = vip_hci.Ds9Window()
    ds9.display(cube[0,0])

psf_norm, flux, fwhm = vip_hci.metrics.normalize_psf(psf,fwhm='fit',full_output=True)


# Stellar photometry of the companion

## Collapse the images for better photometry measurement
cube_derot = vip_hci.preproc.cube_derotate(cube,angs) # Rotate the images to the same north
cube_wl_coll = vip_hci.preproc.cube_collapse(cube_derot,wl_cube=True) # Collapse along the rotation axis - 3D image
cube_coll = vip_hci.preproc.cube_collapse(cube_derot,wl_cube=False) # Collapse along the wavelength axis - 2D image

## Check the collapsed data cubes
if see_collapsed_cube == True:
    ds9 = vip_hci.Ds9Window()
    ds9.display(cube_wl_coll[0],cube_coll) # cube_wl_coll on the left and cube_coll on the right

## Aperture photometry of companions and PSF

### Define photometry
noise_phot = np.zeros_like(wl) #Noise photometry
psf_final_sum = np.zeros_like(wl) #PSF photometry
final_sum = np.zeros_like(wl) #Companion photometry

### Apertures
aper_noise_comp = photutils.CircularAnnulus((145,145),noise_aperture_pos_comp[0],noise_aperture_pos_comp[1])
aper_noise_psf = photutils.CircularAnnulus(psf_pos,noise_aperture_pos_psf[0],noise_aperture_pos_psf[1])

### Aperture photometry
for i in range(0,wl.shape[0]):
    ### Apertures dependent on channel
    aper_comp = photutils.CircularAperture(comp_pos, 1./2*fwhm[i])
    aper_psf = photutils.CircularAperture(psf_pos, 1./2*fwhm[i])
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
                psf_temp = vip_hci.metrics.normalize_psf(psf[i], size=int(13), fwhm=fwhm[i], verbose=False,full_output=True) #3*fwhm
                if i==0:
                    psf_crop = np.zeros((len(psf),psf_temp[0].shape[0],psf_temp[0].shape[1]))
                psf_crop[i] = psf_temp
    cube_negfc = vip_hci.metrics.cube_inject_companions(cube,psf_crop,-angs,flevel=-105,plsc=pxscale,rad_dists=rad_dists,theta=PA) # Remove companion using NEGFC technique

    print("Companion removed")
    print("Computing contrast curve...")
    contrcurve = vip_hci.metrics.contrast_curve(cube_negfc,-angs,psf,4.,pxscale,psf_final_sum,vip_hci.pca.pca,nbranch=n_branches,
              dpi=300, student=False, debug=True ,plot=True, verbose=True, full_output=True, ncomp=ncomp_pca, scale_list=wl)

elif contrast_curves == False:
    print("No contrast curve")
