############################
# Date: 07/08/2019
# Title: Running script for SHIPS for IFS data
# Description: Use this script to run SHIPS for IFS data. In this script you'll find all the necessary parameters to run SHIPS. ONLY SPHERE-DC DATA FOR NOW. VIP and pyKLIP are used.
# VIP version: 0.9.11 (Rainot edit.)
# pyKLIP version: 1.1 NOT IMPLEMENTED YET
# Python version: 3 ONLY
############################

# Set up your parameters

## Define images to analyse
cube_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/ifs_sortframes_dc-IFS_SCIENCE_REDUCED_SPECTRAL_MASTER_CUBE_SORTED-center_im_sorted.fits'
wavelength_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/ifs_sortframes_dc-IFS_SCIENCE_LAMBDA_INFO-lam.fits'
angles_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/ifs_sortframes_dc-IFS_SCIENCE_PARA_ROTATION_CUBE_SORTED-rotnth_sorted.fits'
psf_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/corrected_psf.fits'

# wavelength_filepath = '/home/alan/data/Backup_macbook/SPHERE/IFS/HD93403/ifs_sortframes_dc-IFS_SCIENCE_LAMBDA_INFO-lam.fits'
# cube_filepath = '/home/alan/data/Backup_macbook/SPHERE/IFS/HD93403/ifs_sortframes_dc-IFS_SCIENCE_REDUCED_SPECTRAL_MASTER_CUBE_SORTED-center_im_sorted.fits'
# angles_filepath = '/home/alan/data/Backup_macbook/SPHERE/IFS/HD93403/ifs_sortframes_dc-IFS_SCIENCE_PARA_ROTATION_CUBE_SORTED-rotnth_sorted.fits'
# psf_filepath = '/home/alan/data/Backup_macbook/SPHERE/IFS/HD93403/ifs_sortframes_dc-IFS_SCIENCE_PSF_MASTER_CUBE-median_unsat.fits'

## Photometry
comp_pos = (129,169) # Companion position in pixels from the center of the frame (X,Y)
psf_pos = (33, 33) # PSF position in pixels (X,Y)
radial_dist = 28. # Radial distance of companion in pixels
position_angle = 121.  # Position angle of companion in degrees
noise_aperture_pos_comp = (15,40) # Position in pixels of the circular annulus aperture for noise measurement in the case of the companion
noise_aperture_pos_psf = (12,22) # Position in pixels of the circular annulus aperture for noise measurement in the case of the PSF
size_psf = 31 # What size PSF would you like to use? ODD VALUE ONLY!!

## Computing power
ncores = 4 # Number of cores you are willing to share for the computation

## Do you want to see the image?
see_cube = False # Original cube
see_collapsed_cube = True # Collapsed cube
see_psf_norm = False # Normalised PSF
see_cube_centre = False # Check if the image is centered correctly

## PCA
ncomp_pca = 0 # Number of principal components for PCA
opti_pca = False # Optimise the number of PCA components?
source_pca = (82.,116.) # Source where to optimise the PCA

## SNR maps
snr_maps = False # Would you like to make and save an SNR map to disk?
snr_map_file = '/home/alan/data/Backup_macbook/SPHERE/IFS/QZCardone/SNRmap_VIP11.fits' # Finish the file with .fits

## Detection
adi_frame = False # Would you like to apply ADI on the frame?
adi_plot = False # Would you like to see the resulting plot?
adi_min_scale = -1 # Minimum colour scale for the ADI plot
adi_max_scale = 1 # Maximum colour scale for the ADI plot
detection = False # Would you like the algorithm to detect sources for you? !! WARNING: this is a simple detection !!

detect_sigma = 5 # What sigma limit would you like for the detection?

## Contrast curves
contrast_curves = False # True or False !! computationally intensive !!
n_branches = 1 # Number of branches for contrast curves


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
import random
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
#psf = np.median(psf, axis=1) # Take the median of all PSFs
psf_scaled = np.zeros_like(psf) # The psf will need to be scaled
flevel = np.zeros_like(cube[:,0,0,0]) # Flux level for the companion
flevel = np.array(flevel) # Redefinition - why?

## Get FWHM of images & normalised PSF
psf_med = vip_hci.preproc.cosmetics.cube_crop_frames(psf, size_psf, xy=(32, 32), verbose=True, force=True) # Resize the PSF
psf_norm, maxflux, fwhm = vip_hci.metrics.normalize_psf(psf_med, fwhm='fit',size=None, threshold=None,mask_core=None, model='gauss',imlib='opencv',interpolation='lanczos4',force_odd=True,full_output=True,verbose=False) # maxflux is a dummy variable

# Stellar photometry of the companion

## Collapse the images for better photometry measurement
# cube_derot = vip_hci.preproc.cube_derotate(cube,angs) # Rotate the images to the same north
# cube_wl_coll = vip_hci.preproc.cube_collapse(cube_derot,wl_cube=True) # Collapse along the rotation axis - 3D image
# cube_coll = vip_hci.preproc.cube_collapse(cube_derot,wl_cube=False) # Collapse along the wavelength axis - 2D image
cube_wl_coll = np.zeros_like(cube[:,0,:,:])
for i in range(len(wl)):
        cube_wl_coll[i] = vip_hci.hci_postproc.median_sub(cube[i],-angs,fwhm=fwhm[i],verbose=False) # Rotate & collapse along the rotation axis - 3D image

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
    ### Noise
    phot_noise = photutils.aperture_photometry(cube_wl_coll[i], aper_noise_comp)
    noise_phot[i] = np.array(phot_noise['aperture_sum'])
    ### PSF
    phot_psf = photutils.aperture_photometry(psf[i], aper_psf)
    phot_psf_noise = photutils.aperture_photometry(psf[i], aper_noise_psf)
    psf_bkg_mean = phot_psf_noise['aperture_sum'] / aper_noise_psf.area()
    psf_bkg_sum = psf_bkg_mean * aper_psf.area()
    psf_final_sum[i] = phot_psf['aperture_sum'] - psf_bkg_sum
    ### Companion
    phot = photutils.aperture_photometry(cube_wl_coll[i], aper_comp)
    bkg_mean = (phot_noise['aperture_sum']-phot['aperture_sum']) / (aper_noise_comp.area()-aper_comp.area())
    bkg_sum = bkg_mean * aper_comp.area()
    final_sum[i] = phot['aperture_sum'] - bkg_sum

cube_mask = cube # Make a safe copy of the cube that will be used as the masked cube

# Masking
## Ad

R2 = 16**2 # Radius of the star considered
mask = np.zeros_like(cube[:,0,:,:])
mask_value = np.zeros((39,32,32))
shift = -angs+angs[0] # Angular shift
for h in range(len(angs)):
    x=0
    y=0
    for i in range(38,70): # Iterate through all pixels in X/Y positions
        for j in range(100,132):
            r2 = (i-54)**2+(j-116)**2# Calculate the radius of the pixel from the center pixel
            if r2<=R2: # If the pixels are in the circle of radius R
                mask_value[:,y,x] = cube[:,h,j,i]
            y+=1
        y=0
        x+=1
    x=0
    y=0
    for i in range(34,66): # Iterate through all pixels in X/Y positions
        for j in range(133,165):
            r2 = (i-50)**2+(j-149)**2# Calculate the radius of the pixel from the center pixel
            if r2<=R2: # If the pixels are in the circle of radius R
                cube_mask[:,h,j,i] = mask_value[:,y,x]
            y+=1
        y=0
        x+=1

# # Old version -- smaller
# R2 = 11**2 # Radius of the star considered
# mask = np.zeros_like(cube[:,0,:,:])
# mask_value = np.zeros((39,22,22))
# shift = -angs+angs[0] # Angular shift
# for h in range(len(angs)):
#     x=0
#     y=0
#     for i in range(49,71): # Iterate through all pixels in X/Y positions
#         for j in range(91,113):
#             r2 = (i-60)**2+(j-102)**2# Calculate the radius of the pixel from the center pixel
#             if r2<=R2: # If the pixels are in the circle of radius R
#                 mask_value[:,y,x] = cube[:,h,j,i]
#             y+=1
#         y=0
#         x+=1
#     x=0
#     y=0
#     for i in range(38,60): # Iterate through all pixels in X/Y positions
#         for j in range(138,160):
#             r2 = (i-50)**2+(j-149)**2# Calculate the radius of the pixel from the center pixel
#             if r2<=R2: # If the pixels are in the circle of radius R
#                 cube_mask[:,h,j,i] = mask_value[:,y,x]
#             y+=1
#         y=0
#         x+=1


R2 = 32**2 # Radius of the star considered
mask = np.zeros_like(cube[:,0,:,:])
mask_value = np.zeros((39,64,64))
shift = -angs+angs[0] # Angular shift
for h in range(len(angs)):
    x=0
    y=0
    for i in range(209,273): # Iterate through all pixels in X/Y positions
        for j in range(155,219):
            r2 = (i-241)**2+(j-187)**2# Calculate the radius of the pixel from the center pixel
            if r2<=R2: # If the pixels are in the circle of radius R
                mask_value[:,y,x] = cube[:,h,j,i]
            y+=1
        y=0
        x+=1
    x=0
    y=0
    # for w in range(len(wl)):
    #     mask_value[w] = vip_hci.preproc.frame_rotate(mask_value[w],angle=270)
    for i in range(224,288): # Iterate through all pixels in X/Y positions
        for j in range(76,140):
            r2 = (i-266)**2+(j-108)**2# Calculate the radius of the pixel from the center pixel
            if r2<=R2: # If the pixels are in the circle of radius R
                cube_mask[:,h,j,i] = mask_value[:,y,x]
            y+=1
        y=0
        x+=1

# # # Old version -- smaller
# R2 = 24**2 # Radius of the star considered
# mask = np.zeros_like(cube[:,0,:,:])
# mask_value = np.zeros((39,47,47))
# shift = -angs+angs[0] # Angular shift
# for h in range(len(angs)):
#     x=0
#     y=0
#     for i in range(230,277): # Iterate through all pixels in X/Y positions
#         for j in range(144,191):
#             r2 = (i-253)**2+(j-167)**2# Calculate the radius of the pixel from the center pixel
#             if r2<=R2: # If the pixels are in the circle of radius R
#                 mask_value[:,y,x] = cube[:,h,j,i]
#             y+=1
#         y=0
#         x+=1
#     x=0
#     y=0
#     # for w in range(len(wl)):
#     #     mask_value[w] = vip_hci.preproc.frame_rotate(mask_value[w],angle=270)
#     for i in range(242,289): # Iterate through all pixels in X/Y positions
#         for j in range(80,127):
#             r2 = (i-266)**2+(j-104)**2# Calculate the radius of the pixel from the center pixel
#             if r2<=R2: # If the pixels are in the circle of radius R
#                 cube_mask[:,h,j,i] = mask_value[:,y,x]
#             y+=1
#         y=0
#         x+=1

# contr = vip_hci.metrics.contrcurve.contrast_curve(cube_mask, -angs, psf_norm, fwhm=np.average(fwhm),
#                                           pxscale=0.0074, starphot=psf_final_sum, algo=vip_hci.pca.pca, sigma=5, nbranch=1,
#                                           theta=190, inner_rad=1, wedge=(0, 360), fc_snr=100,
#                                           student=True, transmission=None, smooth=True,
#                                           interp_order=2, plot=True, dpi=100, imlib='opencv',
#                                           debug=True, verbose=True, full_output=False, save_plot=None,
#                                           object_name=None, frame_size=None, figsize=(8, 4), ncomp=1,adimsdi='double',
#                                           scale_list=wl/wl[19])
# #
# # # np.savetxt('/home/alan/Desktop/IFS_dist_arcsec.txt', contr['distance_arcsec'], delimiter='   ')
# # # np.savetxt('/home/alan/Desktop/IFS_mag_contr.txt', contr['mag_contr'], delimiter='   ')
# plt.show()

d = vip_hci.hci_postproc.median_sub(cube_mask[0],-angs,fwhm=fwhm[0],verbose=False)
ds9 = vip_hci.Ds9Window()
ds9.display(cube_mask[0,0],d)
#vip_hci.fits.write_fits("/Users/alan/Desktop/cube_free_IFS.fits",cube_mask)
