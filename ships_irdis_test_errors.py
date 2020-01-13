############################
# Date: 08/08/2019
# Title: Running script for SHIPS for IRDIS data
# Description: Use this script to run SHIPS for IRDIS data. In this script you'll find all the necessary parameters to run SHIPS. ONLY SPHERE-DC DATA FOR NOW. VIP and pyKLIP are used.
# VIP version: 0.9.11 (Rainot edit.)
# pyKLIP version: 1.1 NOT IMPLEMENTED YET
# Python version: 3 ONLY
############################

# Set up your parameters

## Define images to analyse
cube_filepath = '/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/ird_convert_dc-IRD_SCIENCE_REDUCED_MASTER_CUBE-center_im.fits'
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

size_psf = 31
## PCA
ncomp_pca = 1 # Number of principal components for PCA
opti_pca = False # Optimise the number of PCA components?
source = (501,525) # Source where to optimise the PCA

## Spectrum extraction with Simplex Nelder-Mead optimisation
extract_spec = True # Will start the simplex Nelder-Mead optimisation for spectrum extraction
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

### Scaling the PSF for normalisation -- SHOULD I JUST TAKE PSF_NORM INSTEAD?
psf_scaled = np.zeros_like(psf)
for i in range (0,len(psf)):
    psf_scaled[i] = psf[i]/psf_final_sum[i]

# Spectrum extraction with NM
if extract_spec == True:

    ## Define some parameters
    f_guess_pl = max(final_sum_K1) # Flux first guess
    f_range_K1 = np.zeros((len(final_sum_K1),200))
    f_range_K2 = np.zeros((len(final_sum_K2),200))
    for i in range(0,len(final_sum_K1)):
        f_range_K1[i] = np.linspace(0.2*np.abs(final_sum_K1[i]),70 *np.abs(final_sum_K1[i]),200)
        f_range_K2[i] = np.linspace(0.2*np.abs(final_sum_K2[i]),70 *np.abs(final_sum_K2[i]),200)

r_K1 = np.array([])
theta_K1 = np.array([])
f_K1 = np.array([])
r_K2 = np.array([])
theta_K2 = np.array([])
f_K2 = np.array([])
std_r_K1 = np.zeros_like(final_sum_K1)
std_theta_K1 = np.zeros_like(final_sum_K1)
std_flux_K1 = np.zeros_like(final_sum_K1)
std_r_K2 = np.zeros_like(final_sum_K1)
std_theta_K2 = np.zeros_like(final_sum_K1)
std_flux_K2 = np.zeros_like(final_sum_K1)

with open('/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/VIP_simplex_K1.txt') as f:
    header1 = f.readline()
    header2 = f.readline()
    header3 = f.readline()
    header4 = f.readline()
    for line in f:
        line = line.strip()
        columns = line.split()
        r_K1 = np.append(r_K1,float(columns[0]))
        theta_K1 = np.append(theta_K1,float(columns[1]))
        f_K1 = np.append(f_K1,float(columns[2]))
f.close
with open('/Users/alan/Documents/PhD/Data/SPHERE/IRDIS/QZCar/VIP_simplex_K2.txt') as f:
    header1 = f.readline()
    header2 = f.readline()
    header3 = f.readline()
    header4 = f.readline()
    for line in f:
        line = line.strip()
        columns = line.split()
        r_K2 = np.append(r_K2,float(columns[0]))
        theta_K2 = np.append(theta_K2,float(columns[1]))
        f_K2 = np.append(f_K2,float(columns[2]))
f.close

## Error estimation
simplex_options = {'xtol':1e-2, 'maxiter':500, 'maxfev':1000} # Set the simplex options
centy, centx = vip_hci.var.frame_center(cube[0,0])
n_samples=20

for i in range(len(r_K1)-3):
    theta_inj= np.zeros(n_samples)
    r_inj=np.zeros(n_samples)
    for j in range(0,n_samples):
       r_inj[j]=r_K1[i]
       theta_inj[j]=(j+1)*(360./(n_samples+1.))+theta_K1[i]-360.
    print("Companion: ",i)

    spectra_mes_K1 = np.zeros(len(r_inj))
    theta_mes_K1= np.zeros_like(theta_inj)
    r_mes_K1= np.zeros_like(r_inj)
    spectra_mes_K2 = np.zeros(len(r_inj))
    theta_mes_K2= np.zeros_like(theta_inj)
    r_mes_K2= np.zeros_like(r_inj)

    for j in range(0,n_samples):
        posy = r_inj[j] * np.sin(np.deg2rad(theta_inj[j])) + centy
        posx = r_inj[j] * np.cos(np.deg2rad(theta_inj[j])) + centx

        cube_inj=vip_hci.metrics.cube_inject_companions(cube, psf_norm, -angs, flevel=f_K1[i], plsc=0.0074,
                                                       rad_dists=r_inj[j], n_branches=1, theta=theta_inj[j],
                                                       imlib='opencv', interpolation='lanczos4', verbose=False)
        planet_xycoord = np.array([[posx,posy]])
        fake_comp_K1 = vip_hci.negfc.simplex_optim.firstguess(cube_inj[0], -angs, psf_norm[0], ncomp=1, plsc=0.01225,
                                                          fwhm=fwhm[0], annulus_width=3, aperture_radius=3,
                                                          planets_xy_coord=planet_xycoord, cube_ref=None,
                                                          svd_mode='lapack',f_range=f_range_K1[i], simplex=True,
                                                          fmerit='sum',scaling=None, simplex_options=simplex_options,
                                                          collapse='median',verbose=False)

        fake_comp_K2 = vip_hci.negfc.simplex_optim.firstguess(cube_inj[1], -angs, psf_norm[1], ncomp=1, plsc=0.01225,
                                                          fwhm=fwhm[1], annulus_width=3, aperture_radius=3,
                                                          planets_xy_coord=planet_xycoord, cube_ref=None,
                                                          svd_mode='lapack',f_range=f_range_K2[i], simplex=True,
                                                          fmerit='sum',scaling=None, simplex_options=simplex_options,
                                                          collapse='median',verbose=False)

        r_mes_K1[j] = fake_comp_K1[0]
        theta_mes_K1[j] = fake_comp_K1[1]
        spectra_mes_K1[j] = fake_comp_K1[2]
        r_mes_K2[j] = fake_comp_K2[0]
        theta_mes_K2[j] = fake_comp_K2[1]
        spectra_mes_K2[j] = fake_comp_K2[2]

    std_r_K1[i]=np.sqrt(np.sum(abs(r_inj - r_mes_K1)**2)/(n_samples-1))
    std_theta_K1[i]=np.sqrt(np.sum(abs(theta_inj - theta_mes_K1)**2)/(n_samples-1))
    std_flux_K1[i]=np.sqrt(np.sum(abs(spectra_mes_K1 - f_K1[i])**2)/(n_samples-1))
    std_r_K2[i]=np.sqrt(np.sum(abs(r_inj - r_mes_K2)**2)/(n_samples-1))
    std_theta_K2[i]=np.sqrt(np.sum(abs(theta_inj - theta_mes_K2)**2)/(n_samples-1))
    std_flux_K2[i]=np.sqrt(np.sum(abs(spectra_mes_K2 - f_K2[i])**2)/(n_samples-1))
    print("Companion: ",i)
    print("Error r_K1: ", std_r_K1[i])
    print("Error theta_K1: ", std_theta_K1[i])
    print("Error flux_K1: ", std_flux_K1[i])
    print("Error r_K2: ", std_r_K2[i])
    print("Error theta_K2: ", std_theta_K2[i])
    print("Error flux_K2: ", std_flux_K2[i])

np.savetxt("Errors_K1.txt", (std_r_K1,std_theta_K1,std_flux_K1), delimiter='   ') # Saves to file
np.savetxt("Errors_K2.txt", (std_r_K2,std_theta_K2,std_flux_K2), delimiter='   ')
