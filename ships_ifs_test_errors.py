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

## PCA
ncomp_pca = 1 # Number of principal components for PCA
opti_pca = False # Optimise the number of PCA components?
source_pca = (82.,116.) # Source where to optimise the PCA

## Spectrum extraction with Simplex Nelder-Mead optimisation
extract_spec = False # Will start the simplex Nelder-Mead optimisation for spectrum extraction
ann_width = 3 # Annulus width of Simplex
aper_radius = 3 # Aperture Radius of PCA

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
#psf = np.median(psf, axis=1) # Take the median value of all psf observations
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
#cube_derot = vip_hci.preproc.cube_derotate(cube,angs) # Rotate the images to the same north
# cube_wl_coll = vip_hci.preproc.cube_collapse(cube_derot,wl_cube=True) # Collapse along the rotation axis - 3D image
#cube_coll = vip_hci.preproc.cube_collapse(cube_derot,wl_cube=False) # Collapse along the wavelength axis - 2D image

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

# Spectrum extraction with NM

## Define some parameters
comp_xycoord = [(comp_pos[0],comp_pos[1])] # Companion coords
f_guess_pl = max(final_sum) # Flux first guess as the maximum value of the flux
f_range = np.linspace(0.1*f_guess_pl,50*f_guess_pl, 200)
p_in = np.array([radial_dist,PA]) # Regroup companion positions
simplex_options = {'xtol':1e-2, 'maxiter':500, 'maxfev':1000} # Set the simplex options
simplex_guess = np.zeros((39,3)) # Set the simplex variable: r, PA, flux

r = np.array([])
theta = np.array([])
f = np.array([])
std_r = np.zeros_like(final_sum)
std_theta = np.zeros_like(final_sum)
std_flux = np.zeros_like(final_sum)

with open('/Users/alan/Documents/PhD/Data/SPHERE/IFS/QZCardone/VIP_simplex.txt') as f:
    header1 = f.readline()
    for line in f:
        line = line.strip()
        columns = line.split()
        r = np.append(r,float(columns[0]))
        theta = np.append(theta,float(columns[1]))
        f = np.append(f,float(columns[2]))

## Error estimation
simplex_options = {'xtol':1e-2, 'maxiter':500, 'maxfev':1000} # Set the simplex options
centy, centx = vip_hci.var.frame_center(cube[0,0])
n_samples=15

# Load Ad & Ab masked cube
cube_filepath = '/Users/alan/Desktop/cube_free_IFS.fits'
cube_mask = vip_hci.fits.open_fits(cube_filepath)

theta_inj= np.zeros(n_samples)
r_inj=np.zeros(n_samples)
for j in range(0,n_samples):
   r_inj[j]=r[0]
   theta_inj[j]=(j+1)*(360./(n_samples+1.))+theta[0]

f[0] = 2.007886941914275311e+02
spectra_mes = np.zeros(shape=(len(wl),len(r_inj)))
theta_mes= np.zeros_like(theta_inj)
r_mes= np.zeros_like(r_inj)

for j in range(0,n_samples):
   posy = r_inj[j] * np.sin(np.deg2rad(theta_inj[j])) + centy
   posx = r_inj[j] * np.cos(np.deg2rad(theta_inj[j])) + centx
   print('Fake companion #: ', j+1)
   print('r, theta = ',r_inj[j], theta_inj[j])

   cube_inj=vip_hci.metrics.cube_inject_companions(cube_mask, psf_norm, -angs, flevel=f, plsc=0.0074,
                                                   rad_dists=r_inj[j], n_branches=1, theta=theta_inj[j],
                                                   imlib='opencv', interpolation='lanczos4', verbose=False)

   for i in range(0,len(wl)):
       planet_xycoord = np.array([[posx,posy]])
       fake_comp = vip_hci.negfc.simplex_optim.firstguess(cube_inj[i], -angs, psf_norm[i], ncomp=1, plsc=0.0074,
                                                          fwhm=fwhm[i], annulus_width=3, aperture_radius=2,
                                                          planets_xy_coord=planet_xycoord, cube_ref=None,
                                                          svd_mode='lapack',f_range=f_range, simplex=True,
                                                          fmerit='sum',scaling=None, simplex_options=simplex_options,
                                                          collapse='median',verbose=False)
       r_mes[j] = fake_comp[0]
       theta_mes[j] = fake_comp[1]
       spectra_mes[i,j] = fake_comp[2]

std_flux=np.zeros_like(wl)
std_r=np.sqrt(np.sum(abs(r_inj - np.mean(r_inj))**2)/(n_samples-1))
std_theta=np.sqrt(np.sum(abs(theta_inj - np.mean(theta_inj))**2)/(n_samples-1))
std_flux=np.zeros_like(wl)
for i in range(0,len(wl)):
   std_flux[i]=np.sqrt(np.sum(abs(spectra_mes[i,:] - np.mean(spectra_mes[i,:]))**2)/(n_samples-1))

print("Error r_K1: ", std_r)
print("Error theta_K1: ", std_theta)
print("Error flux_K1: ", std_flux)

np.savetxt("Errors_IFS.txt", (std_r,std_theta,std_flux), delimiter='   ') # Saves to file
