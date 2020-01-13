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
R2 = 14**2 # Radius of the star considered
mask = np.zeros_like(cube[:,0,:,:])
mask_value = np.zeros((2,28,28))
shift = -angs+angs[0] # Angular shift
for h in range(len(angs)):
    x=0
    y=0
    for i in range(480,508): # Iterate through all pixels in X/Y positions
        for j in range(557,585):
            r2 = (i-494)**2+(j-571)**2# Calculate the radius of the pixel from the center pixel
            if r2<=R2: # If the pixels are in the circle of radius R
                mask_value[:,y,x] = cube[:,h,j,i]
            y+=1
        y=0
        x+=1
    x=0
    y=0
    for i in range(514,542): # Iterate through all pixels in X/Y positions
        for j in range(557,585):
            r2 = (i-528)**2+(j-571)**2# Calculate the radius of the pixel from the center pixel
            if r2<=R2: # If the pixels are in the circle of radius R
                cube_mask[:,h,j,i] = mask_value[:,y,x]
            y+=1
        y=0
        x+=1

# for i in range(512,542): # Iterate through all pixels in X/Y positions
#     for j in range(560,587):
#         X = random.randint(548,558) # Find random coordinates within a predefined box for noise search
#         Y = random.randint(447,464)
#         mask_value = cube[:,:,Y,X] # Find the noise values at the random coordinates
#         cube_mask[:,:,j,i] = mask_value
#
# ## Ab


R2 = 71**2 # Radius of the star considered
mask = np.zeros_like(cube[:,0,:,:])
mask_value = np.zeros((2,142,142))
shift = -angs+angs[0] # Angular shift
for k in range(0,2):
    if k==1:
        for h in range(len(angs)):
            x=0
            y=0
            for i in range(366,508): # Iterate through all pixels in X/Y positions
                for j in range(505,647):
                    r2 = (i-437)**2+(j-576)**2# Calculate the radius of the pixel from the center pixel
                    if r2<=R2: # If the pixels are in the circle of radius R
                        mask_value[k,y,x] = cube[k,h,j,i]
                    y+=1
                y=0
                x+=1
            x=0
            y=0
            for w in range(0,2):
                mask_value[w] = vip_hci.preproc.frame_rotate(mask_value[w],-90)
            for i in range(402,544): # Iterate through all pixels in X/Y positions
                for j in range(366,508):
                    r2 = (i-472)**2+(j-437)**2# Calculate the radius of the pixel from the center pixel
                    if r2<=R2: # If the pixels are in the circle of radius R
                        cube_mask[k,h,j,i] = mask_value[k,y,x]
                    y+=1
                y=0
                x+=1
    else:
        for h in range(len(angs)):
            x=0
            y=0
            for i in range(366,508): # Iterate through all pixels in X/Y positions
                for j in range(505,647):
                    r2 = (i-437)**2+(j-576)**2# Calculate the radius of the pixel from the center pixel
                    if r2<=R2: # If the pixels are in the circle of radius R
                        mask_value[k,y,x] = cube[k,h,j,i]
                    y+=1
                y=0
                x+=1
            x=0
            y=0
            for w in range(0,2):
                mask_value[w] = vip_hci.preproc.frame_rotate(mask_value[w],-90)
            for i in range(402,544): # Iterate through all pixels in X/Y positions
                for j in range(366,508):
                    r2 = (i-472)**2+(j-437)**2# Calculate the radius of the pixel from the center pixel
                    if r2<=R2: # If the pixels are in the circle of radius R
                        cube_mask[k,h,j,i] = mask_value[k,y,x]
                    y+=1
                y=0
                x+=1

# #Old version
# R2 = 29**2 # Radius of the star considered
# mask = np.zeros_like(cube[:,0,:,:])
# mask_value = np.zeros((2,58,58))
# shift = -angs+angs[0] # Angular shift
# for h in range(len(angs)):
#     x=0
#     y=0
#     for i in range(516,574): # Iterate through all pixels in X/Y positions
#         for j in range(412,470):
#             r2 = (i-545)**2+(j-441)**2# Calculate the radius of the pixel from the center pixel
#             if r2<=R2: # If the pixels are in the circle of radius R
#                 mask_value[:,y,x] = cube[:,h,j,i]
#             y+=1
#         y=0
#         x+=1
#     x=0
#     y=0
#     for i in range(447,505): # Iterate through all pixels in X/Y positions
#         for j in range(410,468):
#             r2 = (i-474)**2+(j-439)**2# Calculate the radius of the pixel from the center pixel
#             if r2<=R2: # If the pixels are in the circle of radius R
#                 cube_mask[:,h,j,i] = mask_value[:,y,x]
#             y+=1
#         y=0
#         x+=1
# for i in range(435,506): # Iterate through all pixels in X/Y positions
#     for j in range(409,467):
#         X = random.randint(428,446) # Find random coordinates within a predefined box for noise search
#         Y = random.randint(536,545)
#         mask_value = cube[:,:,Y,X] # Find the noise values at the random coordinates
#         cube_mask[:,:,j,i] = mask_value
#
## E
R2 = 24**2 # Radius of the star considered
mask = np.zeros_like(cube[:,0,:,:])
mask_value = np.zeros((2,48,48))
shift = -angs+angs[0] # Angular shift
for h in range(len(angs)):
    x=0
    y=0
    for i in range(300,344): # Iterate through all pixels in X/Y positions
        for j in range(394,438):
            r2 = (i-322)**2+(j-416)**2# Calculate the radius of the pixel from the center pixel
            if r2<=R2: # If the pixels are in the circle of radius R
                mask_value[:,y,x] = cube[:,h,j,i]
            y+=1
        y=0
        x+=1
    x=0
    y=0
    for i in range(326,370): # Iterate through all pixels in X/Y positions
        for j in range(356,400):
            r2 = (i-348)**2+(j-378)**2# Calculate the radius of the pixel from the center pixel
            if r2<=R2: # If the pixels are in the circle of radius R
                cube_mask[:,h,j,i] = mask_value[:,y,x]
            y+=1
        y=0
        x+=1

## S?
R2 = 15**2 # Radius of the star considered
mask = np.zeros_like(cube[:,0,:,:])
mask_value = np.zeros((2,30,30))
shift = -angs+angs[0] # Angular shift
for h in range(len(angs)):
    x=0
    y=0
    for i in range(384,414): # Iterate through all pixels in X/Y positions
        for j in range(349,379):
            r2 = (i-399)**2+(j-364)**2# Calculate the radius of the pixel from the center pixel
            if r2<=R2: # If the pixels are in the circle of radius R
                mask_value[:,y,x] = cube[:,h,j,i]
            y+=1
        y=0
        x+=1
    x=0
    y=0
    for i in range(412,442): # Iterate through all pixels in X/Y positions
        for j in range(332,362):
            r2 = (i-427)**2+(j-347)**2# Calculate the radius of the pixel from the center pixel
            if r2<=R2: # If the pixels are in the circle of radius R
                cube_mask[:,h,j,i] = mask_value[:,y,x]
            y+=1
        y=0
        x+=1

## S?
R2 = 25**2 # Radius of the star considered
mask = np.zeros_like(cube[:,0,:,:])
mask_value = np.zeros((2,50,50))
shift = -angs+angs[0] # Angular shift
for h in range(len(angs)):
    x=0
    y=0
    for i in range(331,381): # Iterate through all pixels in X/Y positions
        for j in range(317,367):
            r2 = (i-356)**2+(j-342)**2# Calculate the radius of the pixel from the center pixel
            if r2<=R2: # If the pixels are in the circle of radius R
                mask_value[:,y,x] = cube[:,h,j,i]
            y+=1
        y=0
        x+=1
    x=0
    y=0
    for i in range(375,425): # Iterate through all pixels in X/Y positions
        for j in range(287,337):
            r2 = (i-400)**2+(j-312)**2# Calculate the radius of the pixel from the center pixel
            if r2<=R2: # If the pixels are in the circle of radius R
                cube_mask[:,h,j,i] = mask_value[:,y,x]
            y+=1
        y=0
        x+=1

# contr = vip_hci.metrics.contrcurve.contrast_curve(cube_mask[1], -angs, psf_norm[1], fwhm=fwhm[1],
#                                           pxscale=0.01225, starphot=1033893.75, algo=vip_hci.pca.pca, sigma=5, nbranch=1,
#                                           theta=120, inner_rad=1, wedge=(0, 360), fc_snr=100,
#                                           student=True, transmission=None, smooth=True,
#                                           interp_order=2, plot=True, dpi=100, imlib='opencv',
#                                           debug=True, verbose=True, full_output=False, save_plot=None,
#                                           object_name=None, frame_size=None, figsize=(8, 4), ncomp=1,adimsdi='single',
#                                           scale_list=wl[1]/wl[1])
#
# # np.savetxt('/home/alan/Desktop/IRDIS_dist_arcsec_K2.txt', contr['distance_arcsec'], delimiter='   ')
# # np.savetxt('/home/alan/Desktop/IRDIS_mag_contr_K2.txt', contr['mag_contr'], delimiter='   ')
# plt.show()

#
d = vip_hci.hci_postproc.median_sub(cube_mask[1],-angs,fwhm=fwhm[1],verbose=False)
ds9 = vip_hci.Ds9Window()
ds9.display(cube_mask[0,0],d)
